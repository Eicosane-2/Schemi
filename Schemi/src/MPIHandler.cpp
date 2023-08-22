/*
 * MPIHandler.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "MPIHandler.hpp"

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
		throw exception("Nullptr in parallCellVolume.", errors::systemError);
}

const schemi::surfaceField<schemi::vector>& schemi::MPIHandler::cSdR() const
{
	if (ownerSurfaceDeltaR)
		return *ownerSurfaceDeltaR;
	else
		throw exception("Nullptr in ownerSurfaceDeltaR.", errors::systemError);
}

schemi::MPIHandler::MPIHandler(std::size_t mpi_rank_in,
		std::size_t mpi_size_in) noexcept :
		mpi_rank(mpi_rank_in), mpi_size(mpi_size_in), parallel(mpi_size != 0)
{
}

std::pair<schemi::vector, schemi::scalar> schemi::MPIHandler::correctParallelepipedVector(
		const vector & parallelepiped) const noexcept
{
#ifdef MPI_VERSION
	vector resultParallelepiped(parallelepiped);

	const scalar dx = parallelepiped.v()[0] / mpi_size;

	resultParallelepiped.v_r()[0] = dx;

	const scalar pointX000 = dx * mpi_rank;

	return std::pair<vector, scalar>(resultParallelepiped, pointX000);
#else
	return std::pair<vector, scalar>(parallelepiped, 0);
#endif
}

void schemi::MPIHandler::initialiseBuffersSize(
		[[maybe_unused]] const mesh & mesh)
{
#ifdef MPI_VERSION
	if (!arraysInit)
	{
		tailData_to_send.reset(new mpi_scalar[mesh.tailNumber()]);
		pointData_to_send.reset(new mpi_scalar[mesh.pointNumber()]);

		tailData_to_rec.reset(new mpi_scalar[mesh.tailNumber()]);
		pointData_to_rec.reset(new mpi_scalar[mesh.pointNumber()]);

		tailSurfaceEnd = mesh.tailNumber();

		pointSurfaceStart = mesh.tailNumber() + mesh.innerNumber();

		pointSurfaceEnd = mesh.tailNumber() + mesh.innerNumber()
				+ mesh.pointNumber();

		localCellSize[0] = mesh.cellsSize();
		cellSendBuf.reset(new mpi_scalar[mesh.cellsSize()]);

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
		throw exception("Arrays were already initialized.", errors::MPIError);
#endif
}

void schemi::MPIHandler::initializeParalleMeshData(const mesh & meshOb)
{
	parallCellVolume = std::make_unique<volumeField<scalar>>(meshOb, 0);

	std::vector<std::pair<boundaryConditionType, vector>> ownerSurfaceDeltaR_BC(
			parallCellVolume->boundCond().size());

	for (std::size_t i = 0; i < ownerSurfaceDeltaR_BC.size(); ++i)
	{
		ownerSurfaceDeltaR_BC[i].first = parallCellVolume->boundCond()[i].first;
		ownerSurfaceDeltaR_BC[i].second = vector { 0 };
	}

	ownerSurfaceDeltaR = std::make_unique<surfaceField<vector>>(meshOb, vector {
			0 }, ownerSurfaceDeltaR_BC);

	correctBoundaryConditions(*parallCellVolume);
	correctBoundaryConditions(*ownerSurfaceDeltaR);

	for (std::size_t i = 0; i < meshOb.cellsSize(); ++i)
		parallCellVolume->ref_r()[i] = meshOb.cells()[i].V();

	for (std::size_t i = 0; i < meshOb.surfacesSize(); ++i)
	{
		const std::size_t surfOwner = meshOb.surfaceOwner()[i];

		ownerSurfaceDeltaR->ref_r()[i] = meshOb.surfaces()[i].rC()
				- meshOb.cells()[surfOwner].rC();
	}

	correctBoundaryValues(*parallCellVolume);
	correctBoundaryValues(*ownerSurfaceDeltaR);
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] volumeField<scalar> & corField) const noexcept
{
#ifdef MPI_VERSION
	if (mpi_rank != 0)
	{
		for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
		{
			const std::size_t cellIndex = corField.meshRef().surfaceOwner()[i];

			tailData_to_send.get()[i - tailSurfaceStart] =
					corField.ref()[cellIndex];
		}

		MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
		MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::right_process),
		MPI_COMM_WORLD);
	}

	if (mpi_rank != (mpi_size - 1))
	{
		MPI_Recv(pointData_to_rec.get(), pointSurfaceEnd - pointSurfaceStart,
		MPI_DOUBLE, mpi_rank + 1, static_cast<int>(MPITags::right_process),
		MPI_COMM_WORLD,
		MPI_STATUSES_IGNORE);

		for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
			corField.boundCond_r()[i].second = pointData_to_rec.get()[i
					- pointSurfaceStart];
	}

	if (mpi_rank != (mpi_size - 1))
	{
		for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
		{
			const std::size_t cellIndex = corField.meshRef().surfaceOwner()[i];

			pointData_to_send.get()[i - pointSurfaceStart] =
					corField.ref()[cellIndex];
		}

		MPI_Send(pointData_to_send.get(), pointSurfaceEnd - pointSurfaceStart,
		MPI_DOUBLE, mpi_rank + 1, static_cast<int>(MPITags::left_process),
		MPI_COMM_WORLD);
	}

	if (mpi_rank != 0)
	{
		MPI_Recv(tailData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
		MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::left_process),
		MPI_COMM_WORLD,
		MPI_STATUSES_IGNORE);

		for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
			corField.boundCond_r()[i].second = tailData_to_rec.get()[i
					- tailSurfaceStart];
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] volumeField<vector> & corField) const noexcept
{
#ifdef MPI_VERSION
	for (std::size_t j = 0; j < vector::vsize; ++j)
	{
		if (mpi_rank != 0)
		{
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				tailData_to_send.get()[i - tailSurfaceStart] =
						corField.ref()[cellIndex].v()[j];
			}

			MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::right_process),
			MPI_COMM_WORLD);
		}

		if (mpi_rank != (mpi_size - 1))
		{
			MPI_Recv(pointData_to_rec.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::right_process),
					MPI_COMM_WORLD,
					MPI_STATUSES_IGNORE);

			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						pointData_to_rec.get()[i - pointSurfaceStart];
		}

		if (mpi_rank != (mpi_size - 1))
		{
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				pointData_to_send.get()[i - pointSurfaceStart] =
						corField.ref()[cellIndex].v()[j];
			}

			MPI_Send(pointData_to_send.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::left_process),
					MPI_COMM_WORLD);
		}

		if (mpi_rank != 0)
		{
			MPI_Recv(tailData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::left_process),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						tailData_to_rec.get()[i - tailSurfaceStart];
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
		if (mpi_rank != 0)
		{
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				tailData_to_send.get()[i - tailSurfaceStart] =
						corField.ref()[cellIndex].v()[j];
			}

			MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::right_process),
			MPI_COMM_WORLD);
		}

		if (mpi_rank != (mpi_size - 1))
		{
			MPI_Recv(pointData_to_rec.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::right_process),
					MPI_COMM_WORLD,
					MPI_STATUSES_IGNORE);

			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						pointData_to_rec.get()[i - pointSurfaceStart];
		}

		if (mpi_rank != (mpi_size - 1))
		{
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				pointData_to_send.get()[i - pointSurfaceStart] =
						corField.ref()[cellIndex].v()[j];
			}

			MPI_Send(pointData_to_send.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::left_process),
					MPI_COMM_WORLD);
		}

		if (mpi_rank != 0)
		{
			MPI_Recv(tailData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::left_process),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						tailData_to_rec.get()[i - tailSurfaceStart];
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
		if (mpi_rank != 0)
		{
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				tailData_to_send.get()[i - tailSurfaceStart] =
						corField.ref()[cellIndex].v()[j];
			}

			MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::right_process),
			MPI_COMM_WORLD);
		}

		if (mpi_rank != (mpi_size - 1))
		{
			MPI_Recv(pointData_to_rec.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::right_process),
					MPI_COMM_WORLD,
					MPI_STATUSES_IGNORE);

			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						pointData_to_rec.get()[i - pointSurfaceStart];
		}

		if (mpi_rank != (mpi_size - 1))
		{
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
			{
				const std::size_t cellIndex =
						corField.meshRef().surfaceOwner()[i];

				pointData_to_send.get()[i - pointSurfaceStart] =
						corField.ref()[cellIndex].v()[j];
			}

			MPI_Send(pointData_to_send.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::left_process),
					MPI_COMM_WORLD);
		}

		if (mpi_rank != 0)
		{
			MPI_Recv(tailData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::left_process),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						tailData_to_rec.get()[i - tailSurfaceStart];
		}
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] surfaceField<scalar> & corField) const noexcept
{
#ifdef MPI_VERSION
	if (mpi_rank != 0)
	{
		for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
			tailData_to_send.get()[i - tailSurfaceStart] = corField.ref()[i];

		MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
		MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::right_process),
		MPI_COMM_WORLD);
	}

	if (mpi_rank != (mpi_size - 1))
	{
		MPI_Recv(pointData_to_rec.get(), pointSurfaceEnd - pointSurfaceStart,
		MPI_DOUBLE, mpi_rank + 1, static_cast<int>(MPITags::right_process),
		MPI_COMM_WORLD,
		MPI_STATUSES_IGNORE);

		for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
			corField.boundCond_r()[i].second = pointData_to_rec.get()[i
					- pointSurfaceStart];
	}

	if (mpi_rank != (mpi_size - 1))
	{
		for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
			pointData_to_send.get()[i - pointSurfaceStart] = corField.ref()[i];

		MPI_Send(pointData_to_send.get(), pointSurfaceEnd - pointSurfaceStart,
		MPI_DOUBLE, mpi_rank + 1, static_cast<int>(MPITags::left_process),
		MPI_COMM_WORLD);
	}

	if (mpi_rank != 0)
	{
		MPI_Recv(tailData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
		MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::left_process),
		MPI_COMM_WORLD,
		MPI_STATUSES_IGNORE);

		for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
			corField.boundCond_r()[i].second = tailData_to_rec.get()[i
					- tailSurfaceStart];
	}
#endif
}

void schemi::MPIHandler::correctBoundaryValues(
		[[maybe_unused]] surfaceField<vector> & corField) const noexcept
{
#ifdef MPI_VERSION
	for (std::size_t j = 0; j < vector::vsize; ++j)
	{
		if (mpi_rank != 0)
		{
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				tailData_to_send.get()[i - tailSurfaceStart] =
						corField.ref()[i].v()[j];

			MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::right_process),
			MPI_COMM_WORLD);
		}

		if (mpi_rank != (mpi_size - 1))
		{
			MPI_Recv(pointData_to_rec.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::right_process),
					MPI_COMM_WORLD,
					MPI_STATUSES_IGNORE);

			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						pointData_to_rec.get()[i - pointSurfaceStart];
		}

		if (mpi_rank != (mpi_size - 1))
		{
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				pointData_to_send.get()[i - pointSurfaceStart] =
						corField.ref()[i].v()[j];

			MPI_Send(pointData_to_send.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::left_process),
					MPI_COMM_WORLD);
		}

		if (mpi_rank != 0)
		{
			MPI_Recv(tailData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::left_process),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						tailData_to_rec.get()[i - tailSurfaceStart];
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
		if (mpi_rank != 0)
		{
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				tailData_to_send.get()[i - tailSurfaceStart] =
						corField.ref()[i].v()[j];

			MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::right_process),
			MPI_COMM_WORLD);
		}

		if (mpi_rank != (mpi_size - 1))
		{
			MPI_Recv(pointData_to_rec.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::right_process),
					MPI_COMM_WORLD,
					MPI_STATUSES_IGNORE);

			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						pointData_to_rec.get()[i - pointSurfaceStart];
		}

		if (mpi_rank != (mpi_size - 1))
		{
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				pointData_to_send.get()[i - pointSurfaceStart] =
						corField.ref()[i].v()[j];

			MPI_Send(pointData_to_send.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::left_process),
					MPI_COMM_WORLD);
		}

		if (mpi_rank != 0)
		{
			MPI_Recv(tailData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::left_process),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						tailData_to_rec.get()[i - tailSurfaceStart];
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
		if (mpi_rank != 0)
		{
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				tailData_to_send.get()[i - tailSurfaceStart] =
						corField.ref()[i].v()[j];

			MPI_Send(tailData_to_send.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::right_process),
			MPI_COMM_WORLD);
		}

		if (mpi_rank != (mpi_size - 1))
		{
			MPI_Recv(pointData_to_rec.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::right_process),
					MPI_COMM_WORLD,
					MPI_STATUSES_IGNORE);

			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						pointData_to_rec.get()[i - pointSurfaceStart];
		}

		if (mpi_rank != (mpi_size - 1))
		{
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				pointData_to_send.get()[i - pointSurfaceStart] =
						corField.ref()[i].v()[j];

			MPI_Send(pointData_to_send.get(),
					pointSurfaceEnd - pointSurfaceStart,
					MPI_DOUBLE, mpi_rank + 1,
					static_cast<int>(MPITags::left_process),
					MPI_COMM_WORLD);
		}

		if (mpi_rank != 0)
		{
			MPI_Recv(tailData_to_rec.get(), tailSurfaceEnd - tailSurfaceStart,
			MPI_DOUBLE, mpi_rank - 1, static_cast<int>(MPITags::left_process),
			MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);

			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].second.v_r()[j] =
						tailData_to_rec.get()[i - tailSurfaceStart];
		}
	}
#endif
}

void schemi::MPIHandler::gatherField(const volumeField<scalar> & inField,
		std::valarray<scalar> & retField) noexcept
{
#ifdef MPI_VERSION
	for (std::size_t i = 0; i < inField.size(); ++i)
		cellSendBuf.get()[i] = inField.ref()[i];

	MPI_Gather(cellSendBuf.get(), localCellSize[0], MPI_DOUBLE, gathBuf.get(),
			localCellSize[0], MPI_DOUBLE, root,
			MPI_COMM_WORLD);

	if (mpi_rank == root)
		for (std::size_t i = 0; i < retField.size(); ++i)
			retField[i] = gathBuf.get()[i];
#else
	retField = inField.ref();
#endif
}

bool schemi::MPIHandler::isRoot()
{
	return mpi_rank == root;
}
