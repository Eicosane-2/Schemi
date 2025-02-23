/*
 * MPIHandler.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "MPIHandler.hpp"

#include <fstream>
#include <iostream>

#if defined(MPI_ENABLE) && defined(MPI_DEBUG)
int MPI_Init(int*, char***)
{
	return 1;
}
int MPI_Comm_rank(int, int*)
{
	return 2;
}
int MPI_Comm_size(int, int*)
{
	return 3;
}
int MPI_Barrier(int)
{
	return 4;
}
int MPI_Gather(const std::size_t*, int, int, std::size_t*, int, int, int, int)
{
	return 51;
}
int MPI_Gather(const double*, int, int, double*, int, int, int, int)
{
	return 52;
}
int MPI_Send(const double*, int, int, int, int, int)
{
	return 6;
}
int MPI_Recv(double*, int, int, int, int, int, int)
{
	return 7;
}
int MPI_Bcast(double*, int, int, int, int)
{
	return 81;
}
int MPI_Bcast(std::size_t*, int, int, int, int)
{
	return 82;
}
int MPI_Finalize()
{
	return 9;
}
#endif

std::size_t schemi::MPIHandler::totCellNum() const noexcept
{
	return gathBufSize;
}

const schemi::volumeField<schemi::scalar>& schemi::MPIHandler::Vol() const
{
	if (parallCellVolume)
		return *parallCellVolume;
	else
		[[unlikely]]
		throw exception("Nullptr in parallCellVolume.", errors::systemError);
}

const schemi::surfaceField<schemi::vector>& schemi::MPIHandler::cSdR() const
{
	if (ownerSurfaceDeltaR)
		return *ownerSurfaceDeltaR;
	else
		[[unlikely]]
		throw exception("Nullptr in ownerSurfaceDeltaR.", errors::systemError);
}

schemi::MPIHandler::MPIHandler(std::size_t mpi_rank_in, std::size_t mpi_size_in) :
		envNodes(), mpi_rank(mpi_rank_in), mpi_size(mpi_size_in), parallel(
				mpi_size > 1)
{
#ifdef MPI_VERSION
	std::ifstream mpiRanksFile { "./set/MPIRanks.txt" };
	if (mpiRanksFile.is_open())
		std::cout << "./set/MPIRanks.txt is opened." << std::endl;
	else
		[[unlikely]]
		throw std::ifstream::failure("./set/MPIRanks.txt not found.");

	std::string skipBuffer;

	mpiRanksFile >> skipBuffer >> rankDirMax[0] >> rankDirMax[1]
			>> rankDirMax[2];

	mpiRanksFile.close();

	if (rankDirMax[0] * rankDirMax[1] * rankDirMax[2] != mpi_size)
		throw exception("Wrong total number of nodes in directions.",
				errors::initialisationError);

	const auto XYSlice = rankDirMax[0] * rankDirMax[1];

	rankDir[2] = mpi_rank / XYSlice;

	const auto currentSliceRank = mpi_rank % XYSlice;

	rankDir[1] = currentSliceRank / rankDirMax[0];

	rankDir[0] = currentSliceRank % rankDirMax[0];

	if (rankDir[0] > 0)
		envNodes.tail = mpi_rank - 1;
	if (rankDir[0] < rankDirMax[0] - 1)
		envNodes.point = mpi_rank + 1;
	if (rankDir[2] > 0)
		envNodes.bottom = mpi_rank - XYSlice;
	if (rankDir[1] > 0)
		envNodes.right = mpi_rank - rankDirMax[0];
	if (rankDir[1] < rankDirMax[1] - 1)
		envNodes.left = mpi_rank + rankDirMax[0];
	if (rankDir[2] < rankDirMax[2] - 1)
		envNodes.top = mpi_rank + XYSlice;
#endif
}

std::pair<schemi::vector, schemi::vector> schemi::MPIHandler::correctParallelepipedVector(
		const vector & parallelepiped) const noexcept
{
#ifdef MPI_VERSION
	vector resultParallelepiped(parallelepiped);

	const scalar dx = parallelepiped()[0] / rankDirMax[0];
	const scalar dy = parallelepiped()[1] / rankDirMax[1];
	const scalar dz = parallelepiped()[2] / rankDirMax[2];

	resultParallelepiped.r()[0] = dx;
	resultParallelepiped.r()[1] = dy;
	resultParallelepiped.r()[2] = dz;

	const scalar pointX000 = dx * rankDir[0];
	const scalar pointY000 = dy * rankDir[1];
	const scalar pointZ000 = dz * rankDir[2];

	return std::pair<vector, vector>(resultParallelepiped, { pointX000,
			pointY000, pointZ000 });
#else
	return std::pair<vector, vector>(parallelepiped, { 0.0, 0.0, 0.0 });
#endif
}

void schemi::MPIHandler::initialiseBuffersSize(
		[[maybe_unused]] const mesh & mesh_)
{
#ifdef MPI_VERSION
	if (!arraysInit)
	{
		tailData_to_send.reset(new mpi_scalar[mesh_.tailNumber()]);
		pointData_to_send.reset(new mpi_scalar[mesh_.pointNumber()]);

		tailData_to_rec.reset(new mpi_scalar[mesh_.pointNumber()]);
		pointData_to_rec.reset(new mpi_scalar[mesh_.tailNumber()]);

		tailSurfaceEnd = mesh_.tailNumber();

		pointSurfaceStart = mesh_.tailNumber() + mesh_.innerNumber();

		pointSurfaceEnd = mesh_.tailNumber() + mesh_.innerNumber()
				+ mesh_.pointNumber();

		if (mesh_.taskDimension() == dimensions::task2D)
		{
			rightData_to_send.reset(new mpi_scalar[mesh_.rightNumber()]);
			leftData_to_send.reset(new mpi_scalar[mesh_.leftNumber()]);

			rightData_to_rec.reset(new mpi_scalar[mesh_.leftNumber()]);
			leftData_to_rec.reset(new mpi_scalar[mesh_.rightNumber()]);

			rightSurfaceStart = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber();
			rightSurfaceEnd = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.rightNumber();

			leftSurfaceStart = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.rightNumber();
			leftSurfaceEnd = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.rightNumber()
					+ mesh_.leftNumber();
		}
		else if (mesh_.taskDimension() == dimensions::task3D)
		{
			bottomData_to_send.reset(new mpi_scalar[mesh_.bottomNumber()]);
			rightData_to_send.reset(new mpi_scalar[mesh_.rightNumber()]);
			leftData_to_send.reset(new mpi_scalar[mesh_.leftNumber()]);
			topData_to_send.reset(new mpi_scalar[mesh_.topNumber()]);

			bottomData_to_rec.reset(new mpi_scalar[mesh_.topNumber()]);
			rightData_to_rec.reset(new mpi_scalar[mesh_.leftNumber()]);
			leftData_to_rec.reset(new mpi_scalar[mesh_.rightNumber()]);
			topData_to_rec.reset(new mpi_scalar[mesh_.bottomNumber()]);

			bottomSurfaceStart = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber();
			bottomSurfaceEnd = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.bottomNumber();

			rightSurfaceStart = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.bottomNumber();
			rightSurfaceEnd = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.bottomNumber()
					+ mesh_.rightNumber();

			leftSurfaceStart = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.bottomNumber()
					+ mesh_.rightNumber();
			leftSurfaceEnd = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.bottomNumber()
					+ mesh_.rightNumber() + mesh_.leftNumber();

			topSurfaceStart = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.bottomNumber()
					+ mesh_.rightNumber() + mesh_.leftNumber();
			topSurfaceEnd = mesh_.tailNumber() + mesh_.innerNumber()
					+ mesh_.pointNumber() + mesh_.bottomNumber()
					+ mesh_.rightNumber() + mesh_.leftNumber()
					+ mesh_.topNumber();
		}

		localCellSize[0] = mesh_.cellsSize();
		cellSendBuf.reset(new mpi_scalar[mesh_.cellsSize()]);

		std::unique_ptr<std::size_t> cellsInNodes = nullptr;

		if (mpi_rank == root)
			cellsInNodes.reset(new std::size_t[mpi_size]);

		MPI_Gather(localCellSize, 1, schemi_MPI_SIZE, cellsInNodes.get(), 1,
		schemi_MPI_SIZE, root, MPI_COMM_WORLD);

		if (mpi_rank == root)
		{
			for (std::size_t n = 0; n < mpi_size; ++n)
				gathBufSize += cellsInNodes.get()[n];

			std::cout << "Total size of system: " << gathBufSize << std::endl;

			if (gathBufSize / mpi_size != localCellSize[0])
				throw exception("Noneven number of cells in nodes.",
						errors::MPIError);

			gathBuf.reset(new mpi_scalar[gathBufSize]);
		}

		arraysInit = true;

		std::size_t gathBufSizeArr[1];

		if (mpi_rank == root)
			gathBufSizeArr[0] = gathBufSize;

		MPI_Bcast(gathBufSizeArr, 1, schemi_MPI_SIZE, root,
		MPI_COMM_WORLD);

		gathBufSize = gathBufSizeArr[0];
	}
	else
		[[unlikely]]
		throw exception("Arrays were already initialised.", errors::MPIError);
#endif
}

void schemi::MPIHandler::initialiseParalleMeshData(const mesh & mesh_)
{
#ifdef MPI_VERSION
	if ((mesh_.taskDimension() == dimensions::task1D)
			&& (rankDirMax[1] != 1 || rankDirMax[2] != 1))
		throw exception("For 1D task second and third MPI ranks must be 1.",
				errors::MPIError);
	else if ((mesh_.taskDimension() == dimensions::task2D)
			&& (rankDirMax[2] != 1))
		throw exception("For 2D task third MPI rank must be 1.",
				errors::MPIError);
#endif

	parallCellVolume = std::make_unique<volumeField<scalar>>(mesh_, 0);

	std::vector<std::pair<boundaryConditionType, vector>> ownerSurfaceDeltaR_BC(
			parallCellVolume->boundCond().size());

	for (std::size_t i = 0; i < ownerSurfaceDeltaR_BC.size(); ++i)
	{
		ownerSurfaceDeltaR_BC[i].first = parallCellVolume->boundCond()[i].first;
		ownerSurfaceDeltaR_BC[i].second = vector { 0 };
	}

	ownerSurfaceDeltaR = std::make_unique<surfaceField<vector>>(mesh_, vector {
			0 }, ownerSurfaceDeltaR_BC);

	correctBoundaryConditions(*parallCellVolume);
	correctBoundaryConditions(*ownerSurfaceDeltaR);

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		parallCellVolume->r()[i] = mesh_.cells()[i].V();

	for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
	{
		const std::size_t surfOwner = mesh_.surfaceOwner()[i];

		ownerSurfaceDeltaR->r()[i] = mesh_.surfaces()[i].rC()
				- mesh_.cells()[surfOwner].rC();
	}

	correctBoundaryValues(*parallCellVolume);
	correctBoundaryValues(*ownerSurfaceDeltaR);
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] volumeField<scalar> & corField) const noexcept
{
#ifdef MPI_VERSION
	switch (corField.meshRef().taskDimension())
	{
	case dimensions::task3D:
	{
		/* Send bottom surface data */
		if (!(envNodes.bottom < 0))
		{
			for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				bottomData_to_send.get()[i - bottomSurfaceStart] =
						corField()[cellIndex];
			}

			MPI_Send(bottomData_to_send.get(),
					bottomSurfaceEnd - bottomSurfaceStart,
					MPI_DOUBLE, envNodes.bottom,
					static_cast<int>(MPITags::bottom),
					MPI_COMM_WORLD);
		}
		if (!(envNodes.top < 0))
		{
			MPI_Recv(bottomData_to_rec.get(), topSurfaceEnd - topSurfaceStart,
			MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::bottom),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = bottomData_to_rec.get()[i
						- topSurfaceStart];
		}
		/* Send top surface data */
		if (!(envNodes.top < 0))
		{
			for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				topData_to_send.get()[i - topSurfaceStart] =
						corField()[cellIndex];
			}

			MPI_Send(topData_to_send.get(), topSurfaceEnd - topSurfaceStart,
			MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::top),
			MPI_COMM_WORLD);
		}
		if (!(envNodes.bottom < 0))
		{
			MPI_Recv(topData_to_rec.get(),
					bottomSurfaceEnd - bottomSurfaceStart,
					MPI_DOUBLE, envNodes.bottom, static_cast<int>(MPITags::top),
					MPI_COMM_WORLD,
					MPI_STATUSES_IGNORE);

			for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = topData_to_rec.get()[i
						- bottomSurfaceStart];
		}
	}
		[[fallthrough]];
	case dimensions::task2D:
	{
		/* Send right surface data */
		if (!(envNodes.right < 0))
		{
			for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				rightData_to_send.get()[i - rightSurfaceStart] =
						corField()[cellIndex];
			}

			MPI_Send(rightData_to_send.get(),
					rightSurfaceEnd - rightSurfaceStart,
					MPI_DOUBLE, envNodes.right,
					static_cast<int>(MPITags::right),
					MPI_COMM_WORLD);
		}
		if (!(envNodes.left < 0))
		{
			MPI_Recv(rightData_to_rec.get(), leftSurfaceEnd - leftSurfaceStart,
			MPI_DOUBLE, envNodes.left, static_cast<int>(MPITags::right),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = rightData_to_rec.get()[i
						- leftSurfaceStart];
		}
		/* Send left surface data */
		if (!(envNodes.left < 0))
		{
			for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				leftData_to_send.get()[i - leftSurfaceStart] =
						corField()[cellIndex];
			}

			MPI_Send(leftData_to_send.get(), leftSurfaceEnd - leftSurfaceStart,
			MPI_DOUBLE, envNodes.left, static_cast<int>(MPITags::left),
			MPI_COMM_WORLD);
		}
		if (!(envNodes.right < 0))
		{
			MPI_Recv(leftData_to_rec.get(), rightSurfaceEnd - rightSurfaceStart,
			MPI_DOUBLE, envNodes.right, static_cast<int>(MPITags::left),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = leftData_to_rec.get()[i
						- rightSurfaceStart];
		}
	}
		[[fallthrough]];
	default:
	{
		/* Send tail surface data. */
		if (!(envNodes.tail < 0))
		{
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				tailData_to_send.get()[i - tailSurfaceStart] =
						corField()[cellIndex];
			}

			MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, envNodes.tail, static_cast<int>(MPITags::tail),
			MPI_COMM_WORLD);
		}
		if (!(envNodes.point < 0))
		{
			MPI_Recv(tailData_to_rec.get(), pointSurfaceEnd - pointSurfaceStart,
			MPI_DOUBLE, envNodes.point, static_cast<int>(MPITags::tail),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = tailData_to_rec.get()[i
						- pointSurfaceStart];
		}
		/* Send point surface data. */
		if (!(envNodes.point < 0))
		{
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				pointData_to_send.get()[i - pointSurfaceStart] =
						corField()[cellIndex];
			}

			MPI_Send(pointData_to_send.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, envNodes.point,
					static_cast<int>(MPITags::point),
					MPI_COMM_WORLD);
		}
		if (!(envNodes.tail < 0))
		{
			MPI_Recv(pointData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, envNodes.tail, static_cast<int>(MPITags::point),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = pointData_to_rec.get()[i
						- tailSurfaceStart];
		}
	}
		break;
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] volumeField<vector> & corField) const noexcept
{
#ifdef MPI_VERSION
	for (std::size_t j = 0; j < vector::vsize; ++j)
	{
		switch (corField.meshRef().taskDimension())
		{
		case dimensions::task3D:
		{
			/* Send bottom surface data */
			if (!(envNodes.bottom < 0))
			{
				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					bottomData_to_send.get()[i - bottomSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(bottomData_to_send.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.top < 0))
			{
				MPI_Recv(bottomData_to_rec.get(),
						topSurfaceEnd - topSurfaceStart,
						MPI_DOUBLE, envNodes.top,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							bottomData_to_rec.get()[i - topSurfaceStart];
			}
			/* Send top surface data */
			if (!(envNodes.top < 0))
			{
				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					topData_to_send.get()[i - topSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(topData_to_send.get(), topSurfaceEnd - topSurfaceStart,
				MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::top),
				MPI_COMM_WORLD);
			}
			if (!(envNodes.bottom < 0))
			{
				MPI_Recv(topData_to_rec.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::top),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							topData_to_rec.get()[i - bottomSurfaceStart];
			}
		}
			[[fallthrough]];
		case dimensions::task2D:
		{
			/* Send right surface data */
			if (!(envNodes.right < 0))
			{
				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					rightData_to_send.get()[i - rightSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(rightData_to_send.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.left < 0))
			{
				MPI_Recv(rightData_to_rec.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							rightData_to_rec.get()[i - leftSurfaceStart];
			}
			/* Send left surface data */
			if (!(envNodes.left < 0))
			{
				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					leftData_to_send.get()[i - leftSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(leftData_to_send.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.right < 0))
			{
				MPI_Recv(leftData_to_rec.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							leftData_to_rec.get()[i - rightSurfaceStart];
			}
		}
			[[fallthrough]];
		default:
		{
			/* Send tail surface data. */
			if (!(envNodes.tail < 0))
			{
				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					tailData_to_send.get()[i - tailSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(tailData_to_send.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.point < 0))
			{
				MPI_Recv(tailData_to_rec.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							tailData_to_rec.get()[i - pointSurfaceStart];
			}
			/* Send point surface data. */
			if (!(envNodes.point < 0))
			{
				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					pointData_to_send.get()[i - pointSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(pointData_to_send.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.tail < 0))
			{
				MPI_Recv(pointData_to_rec.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							pointData_to_rec.get()[i - tailSurfaceStart];
			}
		}
			break;
		}
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] volumeField<tensor> & corField) const noexcept
{
#ifdef MPI_VERSION
	for (std::size_t j = 0; j < tensor::vsize; ++j)
	{
		switch (corField.meshRef().taskDimension())
		{
		case dimensions::task3D:
		{
			/* Send bottom surface data */
			if (!(envNodes.bottom < 0))
			{
				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					bottomData_to_send.get()[i - bottomSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(bottomData_to_send.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.top < 0))
			{
				MPI_Recv(bottomData_to_rec.get(),
						topSurfaceEnd - topSurfaceStart,
						MPI_DOUBLE, envNodes.top,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							bottomData_to_rec.get()[i - topSurfaceStart];
			}
			/* Send top surface data */
			if (!(envNodes.top < 0))
			{
				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					topData_to_send.get()[i - topSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(topData_to_send.get(), topSurfaceEnd - topSurfaceStart,
				MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::top),
				MPI_COMM_WORLD);
			}
			if (!(envNodes.bottom < 0))
			{
				MPI_Recv(topData_to_rec.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::top),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							topData_to_rec.get()[i - bottomSurfaceStart];
			}
		}
			[[fallthrough]];
		case dimensions::task2D:
		{
			/* Send right surface data */
			if (!(envNodes.right < 0))
			{
				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					rightData_to_send.get()[i - rightSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(rightData_to_send.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.left < 0))
			{
				MPI_Recv(rightData_to_rec.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							rightData_to_rec.get()[i - leftSurfaceStart];
			}
			/* Send left surface data */
			if (!(envNodes.left < 0))
			{
				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					leftData_to_send.get()[i - leftSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(leftData_to_send.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.right < 0))
			{
				MPI_Recv(leftData_to_rec.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							leftData_to_rec.get()[i - rightSurfaceStart];
			}
		}
			[[fallthrough]];
		default:
		{
			/* Send tail surface data. */
			if (!(envNodes.tail < 0))
			{
				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					tailData_to_send.get()[i - tailSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(tailData_to_send.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.point < 0))
			{
				MPI_Recv(tailData_to_rec.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							tailData_to_rec.get()[i - pointSurfaceStart];
			}
			/* Send point surface data. */
			if (!(envNodes.point < 0))
			{
				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					pointData_to_send.get()[i - pointSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(pointData_to_send.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.tail < 0))
			{
				MPI_Recv(pointData_to_rec.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							pointData_to_rec.get()[i - tailSurfaceStart];
			}
		}
			break;
		}
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] volumeField<tensor3> & corField) const noexcept
{
#ifdef MPI_VERSION
	for (std::size_t j = 0; j < tensor3::vsize; ++j)
	{
		switch (corField.meshRef().taskDimension())
		{
		case dimensions::task3D:
		{
			/* Send bottom surface data */
			if (!(envNodes.bottom < 0))
			{
				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					bottomData_to_send.get()[i - bottomSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(bottomData_to_send.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.top < 0))
			{
				MPI_Recv(bottomData_to_rec.get(),
						topSurfaceEnd - topSurfaceStart,
						MPI_DOUBLE, envNodes.top,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							bottomData_to_rec.get()[i - topSurfaceStart];
			}
			/* Send top surface data */
			if (!(envNodes.top < 0))
			{
				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					topData_to_send.get()[i - topSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(topData_to_send.get(), topSurfaceEnd - topSurfaceStart,
				MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::top),
				MPI_COMM_WORLD);
			}
			if (!(envNodes.bottom < 0))
			{
				MPI_Recv(topData_to_rec.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::top),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							topData_to_rec.get()[i - bottomSurfaceStart];
			}
		}
			[[fallthrough]];
		case dimensions::task2D:
		{
			/* Send right surface data */
			if (!(envNodes.right < 0))
			{
				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					rightData_to_send.get()[i - rightSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(rightData_to_send.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.left < 0))
			{
				MPI_Recv(rightData_to_rec.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							rightData_to_rec.get()[i - leftSurfaceStart];
			}
			/* Send left surface data */
			if (!(envNodes.left < 0))
			{
				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					leftData_to_send.get()[i - leftSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(leftData_to_send.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.right < 0))
			{
				MPI_Recv(leftData_to_rec.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							leftData_to_rec.get()[i - rightSurfaceStart];
			}
		}
			[[fallthrough]];
		default:
		{
			/* Send tail surface data. */
			if (!(envNodes.tail < 0))
			{
				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					tailData_to_send.get()[i - tailSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(tailData_to_send.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.point < 0))
			{
				MPI_Recv(tailData_to_rec.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							tailData_to_rec.get()[i - pointSurfaceStart];
			}
			/* Send point surface data. */
			if (!(envNodes.point < 0))
			{
				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
				{
					const std::size_t cellIndex =
							corField.meshRef().surfaceOwner()[i];

					pointData_to_send.get()[i - pointSurfaceStart] =
							corField()[cellIndex]()[j];
				}

				MPI_Send(pointData_to_send.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.tail < 0))
			{
				MPI_Recv(pointData_to_rec.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							pointData_to_rec.get()[i - tailSurfaceStart];
			}
		}
			break;
		}
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] surfaceField<scalar> & corField) const noexcept
{
#ifdef MPI_VERSION
	switch (corField.meshRef().taskDimension())
	{
	case dimensions::task3D:
	{
		/* Send bottom surface data */
		if (!(envNodes.bottom < 0))
		{
			for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd; ++i)
				bottomData_to_send.get()[i - bottomSurfaceStart] =
						corField()[i];

			MPI_Send(bottomData_to_send.get(),
					bottomSurfaceEnd - bottomSurfaceStart,
					MPI_DOUBLE, envNodes.bottom,
					static_cast<int>(MPITags::bottom),
					MPI_COMM_WORLD);
		}
		if (!(envNodes.top < 0))
		{
			MPI_Recv(bottomData_to_rec.get(), topSurfaceEnd - topSurfaceStart,
			MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::bottom),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = bottomData_to_rec.get()[i
						- topSurfaceStart];
		}
		/* Send top surface data */
		if (!(envNodes.top < 0))
		{
			for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
				topData_to_send.get()[i - topSurfaceStart] = corField()[i];

			MPI_Send(topData_to_send.get(), topSurfaceEnd - topSurfaceStart,
			MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::top),
			MPI_COMM_WORLD);
		}
		if (!(envNodes.bottom < 0))
		{
			MPI_Recv(topData_to_rec.get(),
					bottomSurfaceEnd - bottomSurfaceStart,
					MPI_DOUBLE, envNodes.bottom, static_cast<int>(MPITags::top),
					MPI_COMM_WORLD,
					MPI_STATUSES_IGNORE);

			for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = topData_to_rec.get()[i
						- bottomSurfaceStart];
		}
	}
		[[fallthrough]];
	case dimensions::task2D:
	{
		/* Send right surface data */
		if (!(envNodes.right < 0))
		{
			for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd; ++i)
				rightData_to_send.get()[i - rightSurfaceStart] = corField()[i];

			MPI_Send(rightData_to_send.get(),
					rightSurfaceEnd - rightSurfaceStart,
					MPI_DOUBLE, envNodes.right,
					static_cast<int>(MPITags::right),
					MPI_COMM_WORLD);
		}
		if (!(envNodes.left < 0))
		{
			MPI_Recv(rightData_to_rec.get(), leftSurfaceEnd - leftSurfaceStart,
			MPI_DOUBLE, envNodes.left, static_cast<int>(MPITags::right),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = rightData_to_rec.get()[i
						- leftSurfaceStart];
		}
		/* Send left surface data */
		if (!(envNodes.left < 0))
		{
			for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
				leftData_to_send.get()[i - leftSurfaceStart] = corField()[i];

			MPI_Send(leftData_to_send.get(), leftSurfaceEnd - leftSurfaceStart,
			MPI_DOUBLE, envNodes.left, static_cast<int>(MPITags::left),
			MPI_COMM_WORLD);
		}
		if (!(envNodes.right < 0))
		{
			MPI_Recv(leftData_to_rec.get(), rightSurfaceEnd - rightSurfaceStart,
			MPI_DOUBLE, envNodes.right, static_cast<int>(MPITags::left),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = leftData_to_rec.get()[i
						- rightSurfaceStart];
		}
	}
		[[fallthrough]];
	default:
	{
		/* Send tail surface data. */
		if (!(envNodes.tail < 0))
		{
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				tailData_to_send.get()[i - tailSurfaceStart] = corField()[i];

			MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, envNodes.tail, static_cast<int>(MPITags::tail),
			MPI_COMM_WORLD);
		}
		if (!(envNodes.point < 0))
		{
			MPI_Recv(tailData_to_rec.get(), pointSurfaceEnd - pointSurfaceStart,
			MPI_DOUBLE, envNodes.point, static_cast<int>(MPITags::tail),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = tailData_to_rec.get()[i
						- pointSurfaceStart];
		}
		/* Send point surface data. */
		if (!(envNodes.point < 0))
		{
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				pointData_to_send.get()[i - pointSurfaceStart] = corField()[i];

			MPI_Send(pointData_to_send.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, envNodes.point,
					static_cast<int>(MPITags::point),
					MPI_COMM_WORLD);
		}
		if (!(envNodes.tail < 0))
		{
			MPI_Recv(pointData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, envNodes.tail, static_cast<int>(MPITags::point),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].second = pointData_to_rec.get()[i
						- tailSurfaceStart];
		}
	}
		break;
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] surfaceField<vector> & corField) const noexcept
{
#ifdef MPI_VERSION
	for (std::size_t j = 0; j < vector::vsize; ++j)
	{
		switch (corField.meshRef().taskDimension())
		{
		case dimensions::task3D:
		{
			/* Send bottom surface data */
			if (!(envNodes.bottom < 0))
			{
				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
					bottomData_to_send.get()[i - bottomSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(bottomData_to_send.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.top < 0))
			{
				MPI_Recv(bottomData_to_rec.get(),
						topSurfaceEnd - topSurfaceStart,
						MPI_DOUBLE, envNodes.top,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							bottomData_to_rec.get()[i - topSurfaceStart];
			}
			/* Send top surface data */
			if (!(envNodes.top < 0))
			{
				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
					topData_to_send.get()[i - topSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(topData_to_send.get(), topSurfaceEnd - topSurfaceStart,
				MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::top),
				MPI_COMM_WORLD);
			}
			if (!(envNodes.bottom < 0))
			{
				MPI_Recv(topData_to_rec.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::top),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							topData_to_rec.get()[i - bottomSurfaceStart];
			}
		}
			[[fallthrough]];
		case dimensions::task2D:
		{
			/* Send right surface data */
			if (!(envNodes.right < 0))
			{
				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
					rightData_to_send.get()[i - rightSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(rightData_to_send.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.left < 0))
			{
				MPI_Recv(rightData_to_rec.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							rightData_to_rec.get()[i - leftSurfaceStart];
			}
			/* Send left surface data */
			if (!(envNodes.left < 0))
			{
				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
					leftData_to_send.get()[i - leftSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(leftData_to_send.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.right < 0))
			{
				MPI_Recv(leftData_to_rec.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							leftData_to_rec.get()[i - rightSurfaceStart];
			}
		}
			[[fallthrough]];
		default:
		{
			/* Send tail surface data. */
			if (!(envNodes.tail < 0))
			{
				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
					tailData_to_send.get()[i - tailSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(tailData_to_send.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.point < 0))
			{
				MPI_Recv(tailData_to_rec.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							tailData_to_rec.get()[i - pointSurfaceStart];
			}
			/* Send point surface data. */
			if (!(envNodes.point < 0))
			{
				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
					pointData_to_send.get()[i - pointSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(pointData_to_send.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.tail < 0))
			{
				MPI_Recv(pointData_to_rec.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							pointData_to_rec.get()[i - tailSurfaceStart];
			}
		}
			break;
		}
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] surfaceField<tensor> & corField) const noexcept
{
#ifdef MPI_VERSION
	for (std::size_t j = 0; j < tensor::vsize; ++j)
	{
		switch (corField.meshRef().taskDimension())
		{
		case dimensions::task3D:
		{
			/* Send bottom surface data */
			if (!(envNodes.bottom < 0))
			{
				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
					bottomData_to_send.get()[i - bottomSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(bottomData_to_send.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.top < 0))
			{
				MPI_Recv(bottomData_to_rec.get(),
						topSurfaceEnd - topSurfaceStart,
						MPI_DOUBLE, envNodes.top,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							bottomData_to_rec.get()[i - topSurfaceStart];
			}
			/* Send top surface data */
			if (!(envNodes.top < 0))
			{
				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
					topData_to_send.get()[i - topSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(topData_to_send.get(), topSurfaceEnd - topSurfaceStart,
				MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::top),
				MPI_COMM_WORLD);
			}
			if (!(envNodes.bottom < 0))
			{
				MPI_Recv(topData_to_rec.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::top),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							topData_to_rec.get()[i - bottomSurfaceStart];
			}
		}
			[[fallthrough]];
		case dimensions::task2D:
		{
			/* Send right surface data */
			if (!(envNodes.right < 0))
			{
				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
					rightData_to_send.get()[i - rightSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(rightData_to_send.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.left < 0))
			{
				MPI_Recv(rightData_to_rec.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							rightData_to_rec.get()[i - leftSurfaceStart];
			}
			/* Send left surface data */
			if (!(envNodes.left < 0))
			{
				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
					leftData_to_send.get()[i - leftSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(leftData_to_send.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.right < 0))
			{
				MPI_Recv(leftData_to_rec.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							leftData_to_rec.get()[i - rightSurfaceStart];
			}
		}
			[[fallthrough]];
		default:
		{
			/* Send tail surface data. */
			if (!(envNodes.tail < 0))
			{
				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
					tailData_to_send.get()[i - tailSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(tailData_to_send.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.point < 0))
			{
				MPI_Recv(tailData_to_rec.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							tailData_to_rec.get()[i - pointSurfaceStart];
			}
			/* Send point surface data. */
			if (!(envNodes.point < 0))
			{
				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
					pointData_to_send.get()[i - pointSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(pointData_to_send.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.tail < 0))
			{
				MPI_Recv(pointData_to_rec.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							pointData_to_rec.get()[i - tailSurfaceStart];
			}
		}
			break;
		}
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] surfaceField<tensor3> & corField) const noexcept
{
#ifdef MPI_VERSION
	for (std::size_t j = 0; j < tensor3::vsize; ++j)
	{
		switch (corField.meshRef().taskDimension())
		{
		case dimensions::task3D:
		{
			/* Send bottom surface data */
			if (!(envNodes.bottom < 0))
			{
				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
					bottomData_to_send.get()[i - bottomSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(bottomData_to_send.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.top < 0))
			{
				MPI_Recv(bottomData_to_rec.get(),
						topSurfaceEnd - topSurfaceStart,
						MPI_DOUBLE, envNodes.top,
						static_cast<int>(MPITags::bottom),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							bottomData_to_rec.get()[i - topSurfaceStart];
			}
			/* Send top surface data */
			if (!(envNodes.top < 0))
			{
				for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
					topData_to_send.get()[i - topSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(topData_to_send.get(), topSurfaceEnd - topSurfaceStart,
				MPI_DOUBLE, envNodes.top, static_cast<int>(MPITags::top),
				MPI_COMM_WORLD);
			}
			if (!(envNodes.bottom < 0))
			{
				MPI_Recv(topData_to_rec.get(),
						bottomSurfaceEnd - bottomSurfaceStart,
						MPI_DOUBLE, envNodes.bottom,
						static_cast<int>(MPITags::top),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							topData_to_rec.get()[i - bottomSurfaceStart];
			}
		}
			[[fallthrough]];
		case dimensions::task2D:
		{
			/* Send right surface data */
			if (!(envNodes.right < 0))
			{
				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
					rightData_to_send.get()[i - rightSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(rightData_to_send.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.left < 0))
			{
				MPI_Recv(rightData_to_rec.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::right),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							rightData_to_rec.get()[i - leftSurfaceStart];
			}
			/* Send left surface data */
			if (!(envNodes.left < 0))
			{
				for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
					leftData_to_send.get()[i - leftSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(leftData_to_send.get(),
						leftSurfaceEnd - leftSurfaceStart,
						MPI_DOUBLE, envNodes.left,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.right < 0))
			{
				MPI_Recv(leftData_to_rec.get(),
						rightSurfaceEnd - rightSurfaceStart,
						MPI_DOUBLE, envNodes.right,
						static_cast<int>(MPITags::left),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							leftData_to_rec.get()[i - rightSurfaceStart];
			}
		}
			[[fallthrough]];
		default:
		{
			/* Send tail surface data. */
			if (!(envNodes.tail < 0))
			{
				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
					tailData_to_send.get()[i - tailSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(tailData_to_send.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.point < 0))
			{
				MPI_Recv(tailData_to_rec.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::tail),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
					corField.boundCond_r()[i].second.r()[j] =
							tailData_to_rec.get()[i - pointSurfaceStart];
			}
			/* Send point surface data. */
			if (!(envNodes.point < 0))
			{
				for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd;
						++i)
					pointData_to_send.get()[i - pointSurfaceStart] =
							corField()[i]()[j];

				MPI_Send(pointData_to_send.get(),
						pointSurfaceEnd - pointSurfaceStart,
						MPI_DOUBLE, envNodes.point,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD);
			}
			if (!(envNodes.tail < 0))
			{
				MPI_Recv(pointData_to_rec.get(),
						tailSurfaceEnd - tailSurfaceStart,
						MPI_DOUBLE, envNodes.tail,
						static_cast<int>(MPITags::point),
						MPI_COMM_WORLD,
						MPI_STATUSES_IGNORE);

				for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
					corField.boundCond_r()[i].second.r()[j] =
							pointData_to_rec.get()[i - tailSurfaceStart];
			}
		}
			break;
		}
	}
#endif
}

void schemi::MPIHandler::gatherField(const volumeField<scalar> & inField,
		std::valarray<scalar> & retField) noexcept
{
#ifdef MPI_VERSION
	for (std::size_t i = 0; i < inField.size(); ++i)
		cellSendBuf.get()[i] = inField()[i];

	MPI_Gather(cellSendBuf.get(), localCellSize[0], MPI_DOUBLE, gathBuf.get(),
			localCellSize[0], MPI_DOUBLE, root,
			MPI_COMM_WORLD);

	if (mpi_rank == root)
		for (std::size_t i = 0; i < retField.size(); ++i)
			retField[i] = gathBuf.get()[i];
#else
	retField = inField();
#endif
}

bool schemi::MPIHandler::isRoot()
{
	return mpi_rank == root;
}

const std::array<std::size_t, 3>& schemi::MPIHandler::rankDirection() const noexcept
{
	return rankDir;
}
const std::array<std::size_t, 3>& schemi::MPIHandler::rankDirectionMax() const noexcept
{
	return rankDirMax;
}
