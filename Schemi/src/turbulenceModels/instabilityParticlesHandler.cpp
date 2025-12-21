/*
 * instabilityParticlesHandler.cpp
 *
 *  Created on: 2025/08/25
 *      Author: Maxim Boldyrev
 */

#include "instabilityParticlesHandler.hpp"

#include <iostream>
#include <fstream>
#include <filesystem>

schemi::instabilityParticlesHandler::instabilityParticlesHandler(
		const mesh & meshIn, const MPIHandler & par,
		const volumeField<vector> & uCell, const surfaceField<vector> & uSurf,
		const std::pair<std::size_t, std::string> & readDataPoint) :
		meshRef(meshIn), parallelism(par), boundarySurfaces(0), particlePosition(
				0), cellSurfWeight(0), particlesList(0), particleStatus(0), modelUsed(
				false)
{
	for (std::size_t i = 0; i < meshRef.surfacesSize(); ++i)
		if (meshRef.bndType()[i] != boundaryConditionType::innerSurface)
			boundarySurfaces.push_back(i);

#ifdef MPI_VERSION
	std::ifstream inputFile { "./set/GoncharovModel.txt" };

	if (inputFile.is_open())
		std::cout << "./set/GoncharovModel.txt is opened." << std::endl;
	else
		[[unlikely]]
		throw std::ifstream::failure("./set/GoncharovModel.txt not found.");

	inputFile.precision(ioPrecision);

	std::string modelUsedStringRead;

	inputFile >> modelUsedStringRead;

	modelUsed = onOffMap.at(modelUsedStringRead);

	if (modelUsed)
	{
		std::string skipBuffer;

		inputFile >> skipBuffer;

		inputFile >> listSize;

		if (parallelism.isRoot())
		{
			particlesList.resize(listSize);
			particlePosition.resize(listSize);
			cellSurfWeight.resize(listSize);
			particleStatus.resize(listSize,
					initialisationStatus::notInitialised);
		}

		std::string typeOfInitialisationCriterion;

		inputFile >> skipBuffer;

		inputFile >> typeOfInitialisationCriterion;

		inputFile >> skipBuffer;

		for (std::size_t prt = 0; prt < listSize; ++prt)
			if ((readDataPoint.second == "no")
					|| (readDataPoint.second == "initialisation"))
			{
				vector positionVector_prt;
				std::array<std::size_t, 1> nodeParticleLocated { 0 };

				if (parallelism.isRoot())
				{
					inputFile >> std::get<0>(positionVector_prt.wr())
							>> std::get<1>(positionVector_prt.wr())
							>> std::get<2>(positionVector_prt.wr());
				}

				locateParticleNode(prt, positionVector_prt,
						nodeParticleLocated);

				std::vector<std::size_t> positionData(0);
				std::vector<MPIHandler::mpi_scalar> velocity_arr(0),
						weights_arr(0);
				if (parallelism.isRoot()
						|| (parallelism.mpi_rank
								== std::get<0>(nodeParticleLocated)))
				{
					positionData.resize(2);
					velocity_arr.resize(vector::vsize);
					weights_arr.resize(2);
				}

				if (parallelism.mpi_rank == std::get<0>(nodeParticleLocated))
				{
					positionData[0] = findNearCell(positionVector_prt);
					positionData[1] = findNearSurface(positionVector_prt,
							positionData[0]);

					const std::array<std::size_t, 3> particlePos_prt {
							nodeParticleLocated[0], positionData[0],
							positionData[1] };

					const auto cellSurfWeight_prt = calculateWeights(
							positionVector_prt, particlePos_prt);

					const auto velocity_prt = calculateWeightedValue(uCell,
							uSurf, particlePos_prt, cellSurfWeight_prt);

					for (std::size_t j = 0; j < cellSurfWeight_prt.size(); ++j)
						weights_arr[j] = cellSurfWeight_prt[j];

					for (std::size_t j = 0; j < velocity_prt().size(); ++j)
						velocity_arr[j] = velocity_prt()[j];

					MPI_Send(positionData.data(), 2, schemi_MPI_SIZE,
							parallelism.root,
							static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
							MPI_COMM_WORLD);
					MPI_Send(weights_arr.data(), 2, schemi_MPI_SCALAR,
							parallelism.root,
							static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
							MPI_COMM_WORLD);
					MPI_Send(velocity_arr.data(), vector::vsize,
					schemi_MPI_SCALAR, parallelism.root,
							static_cast<int>(MPIHandler::MPITags::GoncharovVelocity),
							MPI_COMM_WORLD);
				}
				else if (parallelism.isRoot())
				{
					MPI_Recv(positionData.data(), 2, schemi_MPI_SIZE,
							static_cast<int>(std::get<0>(nodeParticleLocated)),
							static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
							MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
					MPI_Recv(weights_arr.data(), 2, schemi_MPI_SCALAR,
							static_cast<int>(std::get<0>(nodeParticleLocated)),
							static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
							MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
					MPI_Recv(velocity_arr.data(), vector::vsize,
					schemi_MPI_SCALAR,
							static_cast<int>(std::get<0>(nodeParticleLocated)),
							static_cast<int>(MPIHandler::MPITags::GoncharovVelocity),
							MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

					std::get<positionType::cell>(particlePosition[prt]) =
							positionData[0];
					std::get<positionType::surface>(particlePosition[prt]) =
							positionData[1];

					cellSurfWeight[prt] = { weights_arr[0], weights_arr[1] };

					std::size_t substance1, substance2;
					int pertType;
					scalar eta_0, lambda, rad;

					inputFile >> substance1 >> substance2 >> pertType >> eta_0
							>> lambda >> rad;

					/*Read turbulent parameters.*/
					scalar CmuIn, C0In, C2In, C3In, Ca1In, Cb1In;
					{
						std::ifstream turbulentParametersFile {
								"./set/turbulentParameters.txt" };
						if (turbulentParametersFile.is_open())
							std::cout
									<< "./set/turbulentParameters.txt is opened."
									<< std::endl;
						else
							[[unlikely]]
							throw std::ifstream::failure(
									"./set/turbulentParameters.txt not found.");

						scalar value;

						turbulentParametersFile >> skipBuffer >> CmuIn;
						for (std::size_t i = 0; i < 12; ++i)
						{
							turbulentParametersFile >> skipBuffer >> value;
						}
						turbulentParametersFile >> skipBuffer >> C0In;
						turbulentParametersFile >> skipBuffer >> value;
						turbulentParametersFile >> skipBuffer >> C2In;
						turbulentParametersFile >> skipBuffer >> C3In;
						turbulentParametersFile >> skipBuffer >> value;
						turbulentParametersFile >> skipBuffer >> Ca1In;
						turbulentParametersFile >> skipBuffer >> Cb1In;

						turbulentParametersFile.close();
					}

					particlesList[prt] = BHRGoncharovTracerModel(
							typeOfInitialisationCriterion, positionVector_prt,
							vector(velocity_arr[0], velocity_arr[1],
									velocity_arr[2]), --substance1,
							--substance2, pertType, eta_0, lambda, rad, CmuIn,
							C0In, C2In, C3In, Ca1In, Cb1In);

					insertCaption(prt);
				}
			}
			else
			{
				vector positionVector_prt;
				std::array<std::size_t, 1> nodeParticleLocated { 0 };

				if (parallelism.isRoot())
				{
					auto readData = dataStream(prt, readDataPoint);

					std::string vectorData1, vectorData2, vectorData3;

					readData >> vectorData1 >> vectorData2 >> vectorData3;
					positionVector_prt = vector(std::stod(vectorData1),
							std::stod(vectorData2), std::stod(vectorData3));
				}

				locateParticleNode(prt, positionVector_prt,
						nodeParticleLocated);

				std::vector<std::size_t> positionData(0);
				std::vector<MPIHandler::mpi_scalar> weights_arr(0);
				if (parallelism.isRoot()
						|| (parallelism.mpi_rank
								== std::get<0>(nodeParticleLocated)))
				{
					positionData.resize(2);
					weights_arr.resize(2);
				}

				if (parallelism.mpi_rank == std::get<0>(nodeParticleLocated))
				{
					positionData[0] = findNearCell(positionVector_prt);
					positionData[1] = findNearSurface(positionVector_prt,
							positionData[0]);

					const std::array<std::size_t, 3> particlePos_prt {
							nodeParticleLocated[0], positionData[0],
							positionData[1] };

					const auto cellSurfWeight_prt = calculateWeights(
							positionVector_prt, particlePos_prt);

					for (std::size_t j = 0; j < cellSurfWeight_prt.size(); ++j)
						weights_arr[j] = cellSurfWeight_prt[j];

					MPI_Send(positionData.data(), 2, schemi_MPI_SIZE,
							parallelism.root,
							static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
							MPI_COMM_WORLD);
					MPI_Send(weights_arr.data(), 2, schemi_MPI_SCALAR,
							parallelism.root,
							static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
							MPI_COMM_WORLD);
				}
				else if (parallelism.isRoot())
				{
					MPI_Recv(positionData.data(), 2, schemi_MPI_SIZE,
							static_cast<int>(std::get<0>(nodeParticleLocated)),
							static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
							MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
					MPI_Recv(weights_arr.data(), 2, schemi_MPI_SCALAR,
							static_cast<int>(std::get<0>(nodeParticleLocated)),
							static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
							MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

					std::get<positionType::cell>(particlePosition[prt]) =
							positionData[0];
					std::get<positionType::surface>(particlePosition[prt]) =
							positionData[1];

					cellSurfWeight[prt] = { weights_arr[0], weights_arr[1] };

					std::size_t substance1, substance2;
					int pertType;
					scalar eta_0, lambda, rad;

					inputFile >> skipBuffer >> skipBuffer >> skipBuffer;

					inputFile >> substance1 >> substance2 >> pertType >> eta_0
							>> lambda >> rad;

					auto readData = dataStream(prt, readDataPoint);

					std::string vectorData1, vectorData2, vectorData3;

					readData >> vectorData1 >> vectorData2 >> vectorData3; // positionVector_prt is already read.

					readData >> vectorData1 >> vectorData2 >> vectorData3;
					const vector inPosition1(std::stod(vectorData1),
							std::stod(vectorData2), std::stod(vectorData3));

					std::array<vector, 4> inVelocity;
					for (std::size_t j = 0; j < inVelocity.size(); ++j)
					{
						readData >> vectorData1 >> vectorData2 >> vectorData3;
						inVelocity[j] = vector(std::stod(vectorData1),
								std::stod(vectorData2), std::stod(vectorData3));
					}

					std::size_t inStep;
					readData >> inStep;

					scalar etaCur, eta1Cur, eta2Cur, rho1Cur, rho2Cur, k0Cur,
							eps0Cur, b0Cur;
					vector a0Cur;
					int curStatus;

					readData >> etaCur >> eta1Cur >> eta2Cur >> rho1Cur
							>> rho2Cur >> k0Cur >> eps0Cur >> b0Cur;
					readData >> vectorData1 >> vectorData2 >> vectorData3;
					a0Cur = vector(std::stod(vectorData1),
							std::stod(vectorData2), std::stod(vectorData3));
					readData >> curStatus;

					/*Read turbulent parameters.*/
					scalar CmuIn, C0In, C2In, C3In, Ca1In, Cb1In;
					{
						std::ifstream turbulentParametersFile {
								"./set/turbulentParameters.txt" };
						if (turbulentParametersFile.is_open())
							std::cout
									<< "./set/turbulentParameters.txt is opened."
									<< std::endl;
						else
							[[unlikely]]
							throw std::ifstream::failure(
									"./set/turbulentParameters.txt not found.");

						scalar value;

						turbulentParametersFile >> skipBuffer >> CmuIn;
						for (std::size_t i = 0; i < 12; ++i)
						{
							turbulentParametersFile >> skipBuffer >> value;
						}
						turbulentParametersFile >> skipBuffer >> C0In;
						turbulentParametersFile >> skipBuffer >> value;
						turbulentParametersFile >> skipBuffer >> C2In;
						turbulentParametersFile >> skipBuffer >> C3In;
						turbulentParametersFile >> skipBuffer >> value;
						turbulentParametersFile >> skipBuffer >> Ca1In;
						turbulentParametersFile >> skipBuffer >> Cb1In;

						turbulentParametersFile.close();
					}

					particlesList[prt] = BHRGoncharovTracerModel(
							typeOfInitialisationCriterion, positionVector_prt,
							inPosition1, inVelocity, inStep, --substance1,
							--substance2, pertType, eta_0, lambda, rad, etaCur,
							eta1Cur, eta2Cur, rho1Cur, rho2Cur, k0Cur, eps0Cur,
							b0Cur, a0Cur, interfaceStatus(curStatus), CmuIn,
							C0In, C2In, C3In, Ca1In, Cb1In);

					clearLines(prt, readDataPoint);
				}
			}
	}

	inputFile.close();
#else
	std::ifstream inputFile { "./set/GoncharovModel.txt" };

	if (inputFile.is_open())
		std::cout << "./set/GoncharovModel.txt is opened." << std::endl;
	else
		[[unlikely]]
		throw std::ifstream::failure("./set/GoncharovModel.txt not found.");

	inputFile.precision(ioPrecision);

	std::string modelUsedStringRead;

	inputFile >> modelUsedStringRead;

	modelUsed = onOffMap.at(modelUsedStringRead);

	if (modelUsed)
	{
		std::string skipBuffer;

		inputFile >> skipBuffer;

		std::size_t numberOfParticles;

		inputFile >> numberOfParticles;

		particlesList.resize(numberOfParticles);
		particlePosition.resize(numberOfParticles);
		cellSurfWeight.resize(numberOfParticles);
		particleStatus.resize(numberOfParticles,
				initialisationStatus::notInitialised);

		std::string typeOfInitialisationCriterion;

		inputFile >> skipBuffer;

		inputFile >> typeOfInitialisationCriterion;

		inputFile >> skipBuffer;

		for (std::size_t prt = 0; prt < particlesList.size(); ++prt)
			if ((readDataPoint.second == "no")
					|| (readDataPoint.second == "initialisation"))
			{
				vector positionVector_prt;

				inputFile >> std::get<0>(positionVector_prt.wr())
						>> std::get<1>(positionVector_prt.wr())
						>> std::get<2>(positionVector_prt.wr());

				if (!checkPosition(positionVector_prt))
					throw exception("Particle out of domain.",
							errors::initialisationError);

				auto & particlePos_prt = particlePosition[prt];

				std::get<positionType::rank>(particlePos_prt) =
						parallelism.mpi_rank;
				std::get<positionType::cell>(particlePos_prt) = findNearCell(
						positionVector_prt);
				std::get<positionType::surface>(particlePos_prt) =
						findNearSurface(positionVector_prt,
								std::get<positionType::cell>(particlePos_prt));

				cellSurfWeight[prt] = calculateWeights(positionVector_prt,
						particlePos_prt);

				const auto velocity_prt = calculateWeightedValue(uCell, uSurf,
						particlePos_prt, cellSurfWeight[prt]);

				std::size_t substance1, substance2;
				int pertType;
				scalar eta_0, lambda, rad;

				inputFile >> substance1 >> substance2 >> pertType >> eta_0
						>> lambda >> rad;

				/*Read turbulent parameters.*/
				scalar CmuIn, C0In, C2In, C3In, Ca1In, Cb1In;
				{
					std::ifstream turbulentParametersFile {
							"./set/turbulentParameters.txt" };
					if (turbulentParametersFile.is_open())
						std::cout << "./set/turbulentParameters.txt is opened."
								<< std::endl;
					else
						[[unlikely]]
						throw std::ifstream::failure(
								"./set/turbulentParameters.txt not found.");

					scalar value;

					turbulentParametersFile >> skipBuffer >> CmuIn;
					for (std::size_t i = 0; i < 12; ++i)
					{
						turbulentParametersFile >> skipBuffer >> value;
					}
					turbulentParametersFile >> skipBuffer >> C0In;
					turbulentParametersFile >> skipBuffer >> value;
					turbulentParametersFile >> skipBuffer >> C2In;
					turbulentParametersFile >> skipBuffer >> C3In;
					turbulentParametersFile >> skipBuffer >> value;
					turbulentParametersFile >> skipBuffer >> Ca1In;
					turbulentParametersFile >> skipBuffer >> Cb1In;

					turbulentParametersFile.close();
				}

				particlesList[prt] = BHRGoncharovTracerModel(
						typeOfInitialisationCriterion, positionVector_prt,
						velocity_prt, --substance1, --substance2, pertType,
						eta_0, lambda, rad, CmuIn, C0In, C2In, C3In, Ca1In,
						Cb1In);

				insertCaption(prt);
			}
			else
			{
				std::size_t substance1, substance2;
				int pertType;
				scalar eta_0, lambda, rad;

				inputFile >> skipBuffer >> skipBuffer >> skipBuffer;

				inputFile >> substance1 >> substance2 >> pertType >> eta_0
						>> lambda >> rad;

				auto readData = dataStream(prt, readDataPoint);

				std::string vectorData1, vectorData2, vectorData3;

				readData >> vectorData1 >> vectorData2 >> vectorData3;
				const vector positionVector_prt(std::stod(vectorData1),
						std::stod(vectorData2), std::stod(vectorData3));

				if (!checkPosition(positionVector_prt))
					throw exception("Particle out of domain.",
							errors::initialisationError);

				auto & particlePos_prt = particlePosition[prt];

				std::get<positionType::rank>(particlePos_prt) =
						parallelism.mpi_rank;
				std::get<positionType::cell>(particlePos_prt) = findNearCell(
						positionVector_prt);
				std::get<positionType::surface>(particlePos_prt) =
						findNearSurface(positionVector_prt,
								std::get<positionType::cell>(particlePos_prt));

				cellSurfWeight[prt] = calculateWeights(positionVector_prt,
						particlePos_prt);

				readData >> vectorData1 >> vectorData2 >> vectorData3;
				const vector inPosition1(std::stod(vectorData1),
						std::stod(vectorData2), std::stod(vectorData3));

				std::array<vector, 4> inVelocity;
				for (std::size_t j = 0; j < inVelocity.size(); ++j)
				{
					readData >> vectorData1 >> vectorData2 >> vectorData3;
					inVelocity[j] = vector(std::stod(vectorData1),
							std::stod(vectorData2), std::stod(vectorData3));
				}

				std::size_t inStep;
				readData >> inStep;

				scalar etaCur, eta1Cur, eta2Cur, rho1Cur, rho2Cur, k0Cur,
						eps0Cur, b0Cur;
				vector a0Cur;
				int curStatus;

				readData >> etaCur >> eta1Cur >> eta2Cur >> rho1Cur >> rho2Cur
						>> k0Cur >> eps0Cur >> b0Cur;
				readData >> vectorData1 >> vectorData2 >> vectorData3;
				a0Cur = vector(std::stod(vectorData1), std::stod(vectorData2),
						std::stod(vectorData3));
				readData >> curStatus;

				/*Read turbulent parameters.*/
				scalar CmuIn, C0In, C2In, C3In, Ca1In, Cb1In;
				{
					std::ifstream turbulentParametersFile {
							"./set/turbulentParameters.txt" };
					if (turbulentParametersFile.is_open())
						std::cout << "./set/turbulentParameters.txt is opened."
								<< std::endl;
					else
						[[unlikely]]
						throw std::ifstream::failure(
								"./set/turbulentParameters.txt not found.");

					scalar value;

					turbulentParametersFile >> skipBuffer >> CmuIn;
					for (std::size_t i = 0; i < 12; ++i)
					{
						turbulentParametersFile >> skipBuffer >> value;
					}
					turbulentParametersFile >> skipBuffer >> C0In;
					turbulentParametersFile >> skipBuffer >> value;
					turbulentParametersFile >> skipBuffer >> C2In;
					turbulentParametersFile >> skipBuffer >> C3In;
					turbulentParametersFile >> skipBuffer >> value;
					turbulentParametersFile >> skipBuffer >> Ca1In;
					turbulentParametersFile >> skipBuffer >> Cb1In;

					turbulentParametersFile.close();
				}

				particlesList[prt] = BHRGoncharovTracerModel(
						typeOfInitialisationCriterion, positionVector_prt,
						inPosition1, inVelocity, inStep, --substance1,
						--substance2, pertType, eta_0, lambda, rad, etaCur,
						eta1Cur, eta2Cur, rho1Cur, rho2Cur, k0Cur, eps0Cur,
						b0Cur, a0Cur, interfaceStatus(curStatus), CmuIn, C0In,
						C2In, C3In, Ca1In, Cb1In);

				clearLines(prt, readDataPoint);
			}
	}

	inputFile.close();
#endif
}

bool schemi::instabilityParticlesHandler::checkPosition(
		const vector & position) const noexcept
{
	// Check, is particle under all boundary surfaces.

	for (std::size_t i = 0; i < boundarySurfaces.size(); ++i)
	{
		const auto surfIndex = boundarySurfaces[i];

		const auto & normale = meshRef.surfaces()[surfIndex].N();
		const auto surfParticleVector = meshRef.surfaces()[surfIndex].rC()
				- position;

		const bool coorientation = (surfParticleVector & normale) >= 0;

		if (!coorientation)
			return false;
	}

	return true;
}

std::size_t schemi::instabilityParticlesHandler::findNearCell(
		const vector & position) noexcept
{
	// Find cell near position for the first time. For constructor and particle changing node.

	std::pair<scalar, std::size_t> near(veryBig, -1);

	for (std::size_t cellIndex = 0; cellIndex < meshRef.cellsSize();
			++cellIndex)
	{
		const auto dr = (position - meshRef.cells()[cellIndex].rC()).mag();

		if (dr < near.first)
			near = { dr, cellIndex };
	}

	return near.second;
}

std::size_t schemi::instabilityParticlesHandler::findNearSurface(
		const vector & position, const std::size_t cellIndex) noexcept
{
	std::pair<scalar, std::size_t> near(veryBig, -1);

	const auto & cellSurfaces = meshRef.surfacesOfCells()[cellIndex];

	for (std::size_t i = 0; i < cellSurfaces.size(); ++i)
	{
		const auto surfIndex = cellSurfaces[i];

		const auto dr = (position - meshRef.surfaces()[surfIndex].rC()).mag();

		if (dr < near.first)
			near = { dr, surfIndex };
	}

	return near.second;
}

std::array<schemi::scalar, 2> schemi::instabilityParticlesHandler::calculateWeights(
		const vector & position, const std::array<std::size_t, 3> & indexes)
{
	const auto nearCellIndex = std::get<1>(indexes);
	const auto nearSurfIndex = std::get<2>(indexes);

	const auto pCell = position - meshRef.cells()[nearCellIndex].rC();
	const auto pSurf = position - meshRef.surfaces()[nearSurfIndex].rC();

	const auto cellLength = pCell.mag();
	const auto surfLength = pSurf.mag();
	const auto sumLength = cellLength + surfLength;

	const auto wei = surfLength / sumLength;

	return
	{	wei, 1 - wei};
}

std::size_t schemi::instabilityParticlesHandler::searchOfNearestCell(
		const vector & newPosition,
		const std::pair<scalar, std::size_t> guessValue,
		std::vector<std::size_t> & checkedCells, std::size_t recursionLevel)
{
	std::vector<std::size_t> newRankNeighb;

	const auto & guessCellCells = meshRef.neighboursOfCells()[guessValue.second];

	for (std::size_t i = 0; i < guessCellCells.size(); ++i)
	{
		const auto cellIndex = guessCellCells[i];

		if (std::find(checkedCells.cbegin(), checkedCells.cend(), cellIndex)
				!= checkedCells.cend())
			continue;

		newRankNeighb.push_back(cellIndex);
	}

	std::pair<scalar, std::size_t> near { veryBig, -1 };

	for (std::size_t i = 0; i < newRankNeighb.size(); ++i)
	{
		const auto cellIndex = newRankNeighb[i];

		const auto dr = (newPosition - meshRef.cells()[cellIndex].rC()).mag();

		if (dr < near.first)
			near = { dr, cellIndex };
	}

	if (near.first > guessValue.first)
		// Nearest cell is the current rank neighbor.
		return guessValue.second;
	else
	{
		if (recursionLevel >= maxRecursionDepthNearestCellSearch)
			throw exception(
					"Search of nearest cell. Recursion level is too high.",
					errors::systemError);
		else
		{
			recursionLevel++;

			checkedCells.insert(checkedCells.cend(), newRankNeighb.cbegin(),
					newRankNeighb.cend());

			searchOfNearestCell(newPosition, near, checkedCells,
					recursionLevel);
		}
	}

	return -1;
}

std::array<std::size_t, 2> schemi::instabilityParticlesHandler::findNewNearPosition(
		const vector & newPosition, const std::size_t oldPositionCellIndex)
{
	std::pair<scalar, std::size_t> near(
			(newPosition - meshRef.cells()[oldPositionCellIndex].rC()).mag(),
			oldPositionCellIndex);

	const auto nearOld(near);

	std::vector<std::size_t> checkedCells;

	// First rank environment cells.
	const auto & firstRankNeis =
			meshRef.neighboursOfCells()[oldPositionCellIndex];

	for (std::size_t i = 0; i < firstRankNeis.size(); ++i)
	{
		const auto cellIndex = firstRankNeis[i];

		const auto dr = (newPosition - meshRef.cells()[cellIndex].rC()).mag();

		if (dr < near.first)
			near = { dr, cellIndex };
	}

	if (near.first > nearOld.first)
	{
		// nearCellIndex does not change, but surface index can.
		const auto newSurfIndex(
				findNearSurface(newPosition, oldPositionCellIndex));

		return
		{	oldPositionCellIndex, newSurfIndex};
	}
	else
	{
		checkedCells.push_back(oldPositionCellIndex);
		checkedCells.insert(checkedCells.cend(), firstRankNeis.cbegin(),
				firstRankNeis.cend());

		const auto newNearCellIndex = searchOfNearestCell(newPosition, near,
				checkedCells, 1);

		const auto newSurfIndex(findNearSurface(newPosition, newNearCellIndex));

		return
		{	newNearCellIndex, newSurfIndex};
	}
}

void schemi::instabilityParticlesHandler::insertCaption(
		const std::size_t particleIndex) const
{
	std::string particleResultName { "./result/GoncharovModelParticle_" };
	particleResultName.append(std::to_string(particleIndex));
	particleResultName.append(".tsv");

	std::ofstream result { particleResultName };

	if (!result.is_open())
		throw std::ofstream::failure(
				"File <<GoncharovModelParticle_>> is not opened.");

	result
			<< "t\tx\ty\tz\teta\tstatus\tk0\teps0\ta0_x\ta0_y\ta0_z\tb\tRe\trelDeltaEta\n";

	result.close();
}

void schemi::instabilityParticlesHandler::clearLines(
		const std::size_t particleIndex,
		const std::pair<std::size_t, std::string> & readDataPoint) const
{
	std::string particleResultName { "./result/GoncharovModelParticle_" };
	particleResultName.append(std::to_string(particleIndex));
	particleResultName.append(".tsv");

	std::ifstream readResult { particleResultName };

	if (!readResult.is_open())
		throw std::ifstream::failure(
				"File <<GoncharovModelParticle_>> is not opened.");

	readResult.precision(ioPrecision);

	std::vector<std::string> linesFromFile;
	std::string buffer;
	for (std::size_t l = 0; l < readDataPoint.first + 2; ++l)
	{
		std::getline(readResult, buffer);
		linesFromFile.push_back(buffer);
	}

	readResult.close();

	std::ofstream result { particleResultName };

	if (!result.is_open())
		throw std::ofstream::failure(
				"File <<GoncharovModelParticle_>> is not opened.");

	result.precision(ioPrecision);

	for (std::size_t l = 0; l < linesFromFile.size(); ++l)
		result << linesFromFile[l] << '\n';

	result.close();
}

std::ifstream schemi::instabilityParticlesHandler::dataStream(
		const std::size_t particleIndex,
		const std::pair<std::size_t, std::string> & readDataPoint) const
{
	const auto nOutputStr = std::to_string(readDataPoint.first);
	const std::size_t length_of_noutput = nOutputStr.size();

	if (lengthOfNumber < length_of_noutput)
		throw exception(
				"Output number length larger than specified number of digits.",
				errors::tooBigOutputNumberError);

	std::string bufOutputN(lengthOfNumber - length_of_noutput, '0');
	bufOutputN.append(nOutputStr);
	std::string fieldDataDirectoryName { "./fieldsOutput/output_" };
	fieldDataDirectoryName.append(bufOutputN);

	if (!std::filesystem::exists(fieldDataDirectoryName))
		throw exception(
				std::string("Directory ") + fieldDataDirectoryName
						+ std::string(" doesn't exits."), errors::systemError);

	fieldDataDirectoryName.append("/");

	const auto fieldDataFileName_prt = fieldDataDirectoryName
			+ std::string("GoncharovModelParticleData_")
			+ std::to_string(particleIndex) + std::string(".dat");

	std::ifstream readData { fieldDataFileName_prt };

	if (!readData.is_open())
		throw std::ifstream::failure(
				"File <<GoncharovModelParticleData_>> is not opened.");

	readData.precision(ioPrecision);

	return readData;
}

void schemi::instabilityParticlesHandler::locateParticleNode(
		[[maybe_unused]] const std::size_t particleIndex,
		[[maybe_unused]] vector & positionVector,
		[[maybe_unused]] std::array<std::size_t, 1> & nodeParticleLocated)
{
#ifdef MPI_VERSION
	std::array<MPIHandler::mpi_scalar, vector::vsize> position_arr {
			std::get<0>(positionVector()), std::get<1>(positionVector()),
			std::get<2>(positionVector()) };

	MPI_Bcast(position_arr.data(), vector::vsize, schemi_MPI_SCALAR,
			parallelism.root,
			MPI_COMM_WORLD);

	if (!parallelism.isRoot())
	{
		std::get<0>(positionVector.wr()) = std::get<0>(position_arr);
		std::get<1>(positionVector.wr()) = std::get<1>(position_arr);
		std::get<2>(positionVector.wr()) = std::get<2>(position_arr);
	}

	const std::array<bool, 1> isPositionMy { checkPosition(positionVector) };

	std::unique_ptr<bool> nodesParticleLocation { nullptr };
	if (parallelism.isRoot())
		nodesParticleLocation.reset(new bool[parallelism.mpi_size]);

	MPI_Gather(isPositionMy.data(), 1, MPI_C_BOOL, nodesParticleLocation.get(),
			1,
			MPI_C_BOOL, parallelism.root, MPI_COMM_WORLD);

	if (parallelism.isRoot())
	{
		std::size_t numberOfEntries { 0 };
		for (std::size_t j = 0; j < parallelism.mpi_size; ++j)
			if (nodesParticleLocation.get()[j])
			{
				std::get<0>(nodeParticleLocated) = j;
				numberOfEntries++;
			}

		if (numberOfEntries > 1)
			throw exception("Particle in more than one domain.",
					errors::initialisationError);
		else if (!numberOfEntries)
			throw exception("Particle is out of all domains.",
					errors::initialisationError);

		std::get<positionType::rank>(particlePosition[particleIndex]) =
				std::get<0>(nodeParticleLocated);
	}

	MPI_Bcast(nodeParticleLocated.data(), 1, schemi_MPI_SIZE, parallelism.root,
	MPI_COMM_WORLD);
#endif
}

void schemi::instabilityParticlesHandler::searchDensity(const std::size_t sub,
		const std::array<std::size_t, 2> & mixSubs,
		const concentrationsPack<cubicCell> & concentrations,
		const std::vector<volumeField<scalar>> & densities,
		const vector & normale, const std::size_t cellIndex, scalar & rhoSub,
		vector & r, std::size_t recursionLevel,
		const boundaryConditionValue & boundVal,
		const std::valarray<scalar> & M) const
{
	const auto x1 = concentrations.v[std::get<0>(mixSubs) + 1].cval()[cellIndex]
			/ concentrations.v[0].cval()[cellIndex];
	const auto x2 = concentrations.v[std::get<1>(mixSubs) + 1].cval()[cellIndex]
			/ concentrations.v[0].cval()[cellIndex];

	if ((x1 + x2) < purityParameter)
		throw exception("There is more than two components in mixture.",
				errors::systemError);

	if (!(((x1 <= 0.5) && (x2 >= 0.5)) || ((x2 <= 0.5) && (x1 >= 0.5))))
		throw exception("Incorrect particle's position.", errors::systemError);

	const auto xSub = concentrations.v[sub + 1].cval()[cellIndex]
			/ concentrations.v[0].cval()[cellIndex];

	if (xSub >= purityParameter)
	{
		rhoSub = densities[sub + 1].cval()[cellIndex];
		r = meshRef.cells()[cellIndex].rC();
		return;
	}

	std::size_t ind1, ind2;

	if (meshRef.neighboursOfCells()[cellIndex].size()
			== meshRef.surfacesOfCells()[cellIndex].size())
	{
		const auto & nearCells = meshRef.neighboursOfCells()[cellIndex];

		std::pair<scalar, std::size_t> maxProj { 0, -1 };
		for (std::size_t i = 0; i < nearCells.size(); ++i)
		{
			const auto outerCellIndex = nearCells[i];

			const auto ccVec = meshRef.cells()[outerCellIndex].rC()
					- meshRef.cells()[cellIndex].rC();

			const auto proj = std::abs(normale & (ccVec / ccVec.mag()));

			if (proj > maxProj.first)
				maxProj = { proj, outerCellIndex };
		}
		ind1 = maxProj.second;

		maxProj = { 0, -1 };
		for (std::size_t i = 0; i < nearCells.size(); ++i)
		{
			const auto outerCellIndex = nearCells[i];

			if (outerCellIndex == ind1)
				continue;

			const auto ccVec = meshRef.cells()[outerCellIndex].rC()
					- meshRef.cells()[cellIndex].rC();

			const auto proj = std::abs(normale & (ccVec / ccVec.mag()));

			if (proj > maxProj.first)
				maxProj = { proj, outerCellIndex };
		}
		ind2 = maxProj.second;

		const auto xSub1 = concentrations.v[sub + 1].cval()[ind1]
				/ concentrations.v[0].cval()[ind1];
		const auto xSub2 = concentrations.v[sub + 1].cval()[ind2]
				/ concentrations.v[0].cval()[ind2];

		if (xSub1 == xSub2)
			throw exception(
					"Can not choose search direction, molar fractions are equal.",
					errors::systemError);

		if (xSub1 >= purityParameter)
		{
			rhoSub = densities[sub + 1].cval()[ind1];
			r = meshRef.cells()[ind1].rC();
			return;
		}
		else if (xSub2 >= purityParameter)
		{
			rhoSub = densities[sub + 1].cval()[ind2];
			r = meshRef.cells()[ind2].rC();
			return;
		}
		else
		{
			if (recursionLevel > maxRecursionDepthRhoSearch)
				throw exception(
						"Search of pure cell. Recursion level is too high.",
						errors::systemError);

			std::size_t nextCell;

			if (xSub2 > xSub1)
				nextCell = ind2;
			else
				nextCell = ind1;

			recursionLevel++;

			searchDensity(sub, mixSubs, concentrations, densities, normale,
					nextCell, rhoSub, r, recursionLevel, boundVal, M);
		}
	}
	else
	{
		const auto & nearSurfaces = meshRef.surfacesOfCells()[cellIndex];

		std::pair<scalar, std::size_t> maxProj { 0, -1 };
		for (std::size_t i = 0; i < nearSurfaces.size(); ++i)
		{
			const auto surfaceIndex = nearSurfaces[i];

			const auto scVec = meshRef.surfaces()[surfaceIndex].rC()
					- meshRef.cells()[cellIndex].rC();

			const auto proj = std::abs(normale & (scVec / scVec.mag()));

			if (proj > maxProj.first)
				maxProj = { proj, surfaceIndex };
		}
		ind1 = maxProj.second;

		maxProj = { 0, -1 };
		for (std::size_t i = 0; i < nearSurfaces.size(); ++i)
		{
			const auto surfaceIndex = nearSurfaces[i];

			if (surfaceIndex == ind1)
				continue;

			const auto scVec = meshRef.surfaces()[surfaceIndex].rC()
					- meshRef.cells()[cellIndex].rC();

			const auto proj = std::abs(normale & (scVec / scVec.mag()));

			if (proj > maxProj.first)
				maxProj = { proj, surfaceIndex };
		}
		ind2 = maxProj.second;

		const auto & concentrationBoundary =
				concentrations.v[sub + 1].boundCond();

		if ((concentrationBoundary[ind1].first
				== boundaryConditionType::innerSurface)
				&& (concentrationBoundary[ind2].first
						== boundaryConditionType::innerSurface))
		{
			const auto cellOwner1 = meshRef.surfaceOwner()[ind1];
			const auto cellNeighb1 = meshRef.surfaceNeighbour()[ind1];

			std::size_t outerCell1;
			if (cellOwner1 == cellIndex)
				outerCell1 = cellNeighb1;
			else if (cellNeighb1 == cellIndex)
				outerCell1 = cellOwner1;
			else
				[[unlikely]]
				throw exception("Cell is not a owner or neighbour.",
						errors::systemError);

			const auto xSub1 = concentrations.v[sub + 1].cval()[outerCell1]
					/ concentrations.v[0].cval()[outerCell1];

			const auto cellOwner2 = meshRef.surfaceOwner()[ind2];
			const auto cellNeighb2 = meshRef.surfaceNeighbour()[ind2];

			std::size_t outerCell2;
			if (cellOwner2 == cellIndex)
				outerCell2 = cellNeighb2;
			else if (cellNeighb2 == cellIndex)
				outerCell2 = cellOwner2;
			else
				[[unlikely]]
				throw exception("Cell is not a owner or neighbour.",
						errors::systemError);

			const auto xSub2 = concentrations.v[sub + 1].cval()[outerCell2]
					/ concentrations.v[0].cval()[outerCell2];

			if (xSub1 == xSub2)
				throw exception(
						"Can not choose search direction, molar fractions are equal.",
						errors::systemError);

			if (xSub1 >= purityParameter)
			{
				rhoSub = densities[sub + 1].cval()[outerCell1];
				r = meshRef.cells()[outerCell1].rC();
				return;
			}
			else if (xSub2 >= purityParameter)
			{
				rhoSub = densities[sub + 1].cval()[outerCell2];
				r = meshRef.cells()[outerCell2].rC();
				return;
			}
			else
			{
				if (recursionLevel > maxRecursionDepthRhoSearch)
					[[unlikely]]
					throw exception(
							"Search of pure cell. Recursion level is too high.",
							errors::systemError);

				std::size_t nextCell;

				if (xSub2 > xSub1)
					nextCell = outerCell2;
				else
					nextCell = outerCell1;

				recursionLevel++;

				searchDensity(sub, mixSubs, concentrations, densities, normale,
						nextCell, rhoSub, r, recursionLevel, boundVal, M);
			}
		}
		else if (concentrationBoundary[ind1].first
				== boundaryConditionType::innerSurface
				&& ((concentrationBoundary[ind2].first
						== boundaryConditionType::calculatedParallelBoundary)
						|| (concentrationBoundary[ind2].first
								== boundaryConditionType::fixedValueCell)))
		{
			const auto cellOwner1 = meshRef.surfaceOwner()[ind1];
			const auto cellNeighb1 = meshRef.surfaceNeighbour()[ind1];

			std::size_t outerCell1;
			if (cellOwner1 == cellIndex)
				outerCell1 = cellNeighb1;
			else if (cellNeighb1 == cellIndex)
				outerCell1 = cellOwner1;
			else
				[[unlikely]]
				throw exception("Cell is not a owner or neighbour.",
						errors::systemError);

			const auto xSub1 = concentrations.v[sub + 1].cval()[outerCell1]
					/ concentrations.v[0].cval()[outerCell1];

			const scalar boundaryConc = boundVal.boundaryConditionValueCell(
					concentrations.v[sub + 1].cval()[cellIndex],
					concentrations.v[sub + 1].boundCond()[ind2], cellIndex,
					ind2);

			scalar sumBoundConc(0);
			for (std::size_t comp = 1; comp < concentrations.v.size(); ++comp)
			{
				sumBoundConc += boundVal.boundaryConditionValueCell(
						concentrations.v[comp].cval()[cellIndex],
						concentrations.v[comp].boundCond()[ind2], cellIndex,
						ind2);
			}

			const auto xSub2 = boundaryConc / sumBoundConc;

			if (xSub1 == xSub2)
				throw exception(
						"Can not choose search direction, molar fractions are equal.",
						errors::systemError);

			if (xSub1 >= purityParameter)
			{
				rhoSub = densities[sub + 1].cval()[outerCell1];
				r = meshRef.cells()[outerCell1].rC();
				return;
			}
			else if (xSub2 >= purityParameter)
			{
				rhoSub = boundaryConc * M[sub];
				r = meshRef.cells()[cellIndex].rC()
						+ (meshRef.surfaces()[ind2].rC()
								- meshRef.cells()[cellIndex].rC()) * 2;
				return;
			}
			else
			{
				if (recursionLevel > maxRecursionDepthRhoSearch)
					[[unlikely]]
					throw exception(
							"Search of pure cell. Recursion level is too high.",
							errors::systemError);

				std::size_t nextCell;

				if (xSub2 > xSub1)
				{
					const scalar deltaX = xSub2 - x2;
					const scalar deltaR = (meshRef.surfaces()[ind2].rC()
							- meshRef.cells()[cellIndex].rC()).mag() * 2;

					const auto dXdR = deltaX / deltaR;

					if (dXdR <= 0)
						[[unlikely]]
						throw exception(
								"Don't know how to extrapolate density.",
								errors::systemError);

					const auto deltaRExtrap = (1 - xSub2) / dXdR;

					const auto deltaRho = boundaryConc * M[sub]
							- densities[sub + 1].cval()[cellIndex];

					const auto extrapRho = deltaRho / deltaR * deltaRExtrap;

					const vector extrapR = (meshRef.surfaces()[ind2].rC()
							- meshRef.cells()[cellIndex].rC()) / (0.5 * deltaR)
							* deltaRExtrap + meshRef.cells()[cellIndex].rC();

					rhoSub = extrapRho;
					r = extrapR;

					return;
				}
				else
					nextCell = outerCell1;

				recursionLevel++;

				searchDensity(sub, mixSubs, concentrations, densities, normale,
						nextCell, rhoSub, r, recursionLevel, boundVal, M);
			}
		}
		else if (((concentrationBoundary[ind1].first
				== boundaryConditionType::calculatedParallelBoundary)
				|| (concentrationBoundary[ind1].first
						== boundaryConditionType::fixedValueCell))
				&& concentrationBoundary[ind2].first
						== boundaryConditionType::innerSurface)
		{
			const scalar boundaryConc = boundVal.boundaryConditionValueCell(
					concentrations.v[sub + 1].cval()[cellIndex],
					concentrations.v[sub + 1].boundCond()[ind1], cellIndex,
					ind1);

			scalar sumBoundConc(0);
			for (std::size_t comp = 1; comp < concentrations.v.size(); ++comp)
			{
				sumBoundConc += boundVal.boundaryConditionValueCell(
						concentrations.v[comp].cval()[cellIndex],
						concentrations.v[comp].boundCond()[ind1], cellIndex,
						ind1);
			}

			const auto xSub1 = boundaryConc / sumBoundConc;

			const auto cellOwner2 = meshRef.surfaceOwner()[ind2];
			const auto cellNeighb2 = meshRef.surfaceNeighbour()[ind2];

			std::size_t outerCell2;
			if (cellOwner2 == cellIndex)
				outerCell2 = cellNeighb2;
			else if (cellNeighb2 == cellIndex)
				outerCell2 = cellOwner2;
			else
				[[unlikely]]
				throw exception("Cell is not a owner or neighbour.",
						errors::systemError);

			const auto xSub2 = concentrations.v[sub + 1].cval()[outerCell2]
					/ concentrations.v[0].cval()[outerCell2];

			if (xSub1 == xSub2)
				throw exception(
						"Can not choose search direction, molar fractions are equal.",
						errors::systemError);

			if (xSub1 >= purityParameter)
			{
				rhoSub = boundaryConc * M[sub];
				r = meshRef.cells()[cellIndex].rC()
						+ (meshRef.surfaces()[ind1].rC()
								- meshRef.cells()[cellIndex].rC()) * 2;
				return;
			}
			else if (xSub2 >= purityParameter)
			{
				rhoSub = densities[sub + 1].cval()[outerCell2];
				r = meshRef.cells()[outerCell2].rC();
				return;
			}
			else
			{
				if (recursionLevel > maxRecursionDepthRhoSearch)
					[[unlikely]]
					throw exception(
							"Search of pure cell. Recursion level is too high.",
							errors::systemError);

				std::size_t nextCell;

				if (xSub1 > xSub2)
				{
					const scalar deltaX = xSub1 - x1;
					const scalar deltaR = (meshRef.surfaces()[ind1].rC()
							- meshRef.cells()[cellIndex].rC()).mag() * 2;

					const auto dXdR = deltaX / deltaR;

					if (dXdR <= 0)
						[[unlikely]]
						throw exception(
								"Don't know how to extrapolate density.",
								errors::systemError);

					const auto deltaRExtrap = (1 - xSub1) / dXdR;

					const auto deltaRho = boundaryConc * M[sub]
							- densities[sub + 1].cval()[cellIndex];

					const auto extrapRho = deltaRho / deltaR * deltaRExtrap;

					const vector extrapR = (meshRef.surfaces()[ind1].rC()
							- meshRef.cells()[cellIndex].rC()) / (0.5 * deltaR)
							* deltaRExtrap + meshRef.cells()[cellIndex].rC();

					rhoSub = extrapRho;
					r = extrapR;

					return;
				}
				else
					nextCell = outerCell2;

				recursionLevel++;

				searchDensity(sub, mixSubs, concentrations, densities, normale,
						nextCell, rhoSub, r, recursionLevel, boundVal, M);
			}
		}
		else
			[[unlikely]]
			throw exception(
					"Could not choose search variant of pure density cell.",
					errors::systemError);
	}
}

schemi::volumeField<schemi::scalar> schemi::instabilityParticlesHandler::formProfile(
		const std::array<std::size_t, 2> & substancesIndex,
		const vector & tracerPosition, const scalar radiusOfInfluence,
		const concentrationsPack<cubicCell> & concentrations,
		const boundaryConditionValue & bnc) const
{
	volumeField<scalar> turbProfileFunction { meshRef };

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
	{
		const scalar x1_i =
				concentrations.v[std::get<0>(substancesIndex) + 1].cval()[i]
						/ concentrations.v[0].cval()[i];
		const scalar x2_i =
				concentrations.v[std::get<1>(substancesIndex) + 1].cval()[i]
						/ concentrations.v[0].cval()[i];

		if (x1_i * x2_i < (1 - purityParameter))
		{
			turbProfileFunction.val()[i] = 0;
			continue;
		}

		const auto & cellR = meshRef.cells()[i].rC();
		if ((tracerPosition - cellR).mag() > radiusOfInfluence)
		{
			turbProfileFunction.val()[i] = 0;
			continue;
		}

		scalar xMax;
		std::size_t xMaxInd;
		if (x1_i >= x2_i)
		{
			xMax = x1_i;
			xMaxInd = std::get<0>(substancesIndex);
		}
		else
		{
			xMax = x2_i;
			xMaxInd = std::get<1>(substancesIndex);
		}

		const auto cellValue = GoncharovTracerModel::C01 * (1 - xMax);

		const auto nearCellSize = meshRef.neighboursOfCells()[i].size();
		const auto cellSurfacesSize = meshRef.surfacesOfCells()[i].size();

		if (nearCellSize != cellSurfacesSize)
		{
			std::valarray<scalar> nearCellValues(0.0, nearCellSize);
			std::valarray<scalar> nearSurfValues(0.0, cellSurfacesSize);

			for (std::size_t j = 0; j < nearCellSize; ++j)
			{
				const auto cellIndex = meshRef.neighboursOfCells()[i][j];

				nearCellValues[j] =
						std::abs(
								xMax
										- concentrations.v[xMaxInd + 1].cval()[cellIndex]
												/ concentrations.v[0].cval()[cellIndex]);
			}

			for (std::size_t l = 0; l < cellSurfacesSize; ++l)
			{
				const auto surfaceIndex = meshRef.surfacesOfCells()[i][l];

				if (concentrations.v[xMaxInd + 1].boundCond()[surfaceIndex].first
						!= boundaryConditionType::innerSurface)
				{
					const auto concMax =
							bnc.boundaryConditionValueCell(
									concentrations.v[xMaxInd + 1].cval()[i],
									concentrations.v[xMaxInd + 1].boundCond()[surfaceIndex],
									i, surfaceIndex);

					scalar sumConc(0);
					for (std::size_t v = 1; v < concentrations.v.size(); ++v)
					{
						sumConc += bnc.boundaryConditionValueCell(
								concentrations.v[v].cval()[i],
								concentrations.v[v].boundCond()[surfaceIndex],
								i, surfaceIndex);
					}

					nearSurfValues[l] = std::abs(xMax - concMax / sumConc);
				}
			}

			turbProfileFunction.val()[i] = std::min(
					cellValue
							+ GoncharovTracerModel::C02
									* (nearCellValues.sum()
											+ nearSurfValues.sum()), 1.0);
		}
		else
		{
			std::valarray<scalar> nearCellValues(0.0, nearCellSize);

			for (std::size_t j = 0; j < nearCellSize; ++j)
			{
				const auto cellIndex = meshRef.neighboursOfCells()[i][j];

				nearCellValues[j] =
						std::abs(
								xMax
										- concentrations.v[xMaxInd + 1].cval()[cellIndex]
												/ concentrations.v[0].cval()[cellIndex]);
			}

			turbProfileFunction.val()[i] = std::min(
					cellValue
							+ GoncharovTracerModel::C02 * nearCellValues.sum(),
					1.0);
		}
	}

	return turbProfileFunction;
}

void schemi::instabilityParticlesHandler::timeIntegration(
		const volumeField<vector> & gradRhoCell,
		const surfaceField<vector> & gradRhoSurf,
		const volumeField<vector> & uCell, const surfaceField<vector> & uSurf,
		const concentrationsPack<cubicCell> & concentrations,
		const std::vector<volumeField<scalar>> & densities,
		const boundaryConditionValue & boundVal,
		const std::valarray<scalar> & M, const scalar timestep,
		const volumeField<vector> & gradP, const volumeField<scalar> & divU,
		const volumeField<tensor> & gradU)
{
#ifdef MPI_VERSION
	for (std::size_t prt = 0; prt < listSize; ++prt)
	{
		std::vector<int> interfStatus_arr(0);
		std::vector<std::size_t> partLocOld_arr(0), substances_arr(0);
		std::vector<MPIHandler::mpi_scalar> weightsOld_arr(0), velocity_arr(0),
				gradRho_arr(0), rho_arr(0), gradP_arr(0), divU_arr(0),
				gradU_arr(0);

		std::array<std::size_t, 1> nodeParticleLocatedOld { 0 };

		if (parallelism.isRoot())
			std::get<0>(nodeParticleLocatedOld) = std::get<positionType::rank>(
					particlePosition[prt]);

		MPI_Bcast(nodeParticleLocatedOld.data(), 1, schemi_MPI_SIZE,
				parallelism.root,
				MPI_COMM_WORLD);

		if (parallelism.isRoot()
				|| (parallelism.mpi_rank == std::get<0>(nodeParticleLocatedOld)))
		{
			interfStatus_arr.resize(1);
			partLocOld_arr.resize(2);
			substances_arr.resize(2);
			weightsOld_arr.resize(2);
			velocity_arr.resize(vector::vsize);
			gradRho_arr.resize(vector::vsize);
			rho_arr.resize(3);
			gradP_arr.resize(vector::vsize);
			divU_arr.resize(1);
			gradU_arr.resize(tensor::vsize);
		}

		if (parallelism.isRoot())
		{
			interfStatus_arr[0] = static_cast<int>(particleStatus[prt]);
			for (std::size_t j = 0; j < 2; ++j)
			{
				partLocOld_arr[j] = particlePosition[prt][j + 1];
				substances_arr[j] = particlesList[prt].substances()[j];
				weightsOld_arr[j] = cellSurfWeight[prt][j];
			}

			MPI_Send(interfStatus_arr.data(), 1, MPI_INT,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovInterfaceStatus),
					MPI_COMM_WORLD);
			MPI_Send(partLocOld_arr.data(), 2, schemi_MPI_SIZE,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
					MPI_COMM_WORLD);
			MPI_Send(substances_arr.data(), 2, schemi_MPI_SIZE,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovSubstances),
					MPI_COMM_WORLD);
			MPI_Send(weightsOld_arr.data(), 2, schemi_MPI_SCALAR,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
					MPI_COMM_WORLD);
		}
		else if (parallelism.mpi_rank == std::get<0>(nodeParticleLocatedOld))
		{
			MPI_Recv(interfStatus_arr.data(), 1, MPI_INT, parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovInterfaceStatus),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(partLocOld_arr.data(), 2, schemi_MPI_SIZE,
					parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(substances_arr.data(), 2, schemi_MPI_SIZE,
					parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovSubstances),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(weightsOld_arr.data(), 2, schemi_MPI_SCALAR,
					parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			const initialisationStatus status { interfStatus_arr[0] };
			const std::array<std::size_t, 3> particlePosition_prt {
					nodeParticleLocatedOld[0], partLocOld_arr[0],
					partLocOld_arr[1] };
			const std::array<std::size_t, 2> particleSubstances_prt {
					substances_arr[0], substances_arr[1] };
			const std::array<scalar, 2> cellSurfWeight_prt { weightsOld_arr[0],
					weightsOld_arr[1] };

			const auto particleVelocity = calculateWeightedValue(uCell, uSurf,
					particlePosition_prt, cellSurfWeight_prt);
			const auto particleDensityGrad = calculateWeightedValue(gradRhoCell,
					gradRhoSurf, particlePosition_prt, cellSurfWeight_prt);

			const auto normale = particleDensityGrad
					/ -particleDensityGrad.mag();

			scalar rho1 { 0 }, rho2 { 0 };
			vector r1 { 0, 0, 0 }, r2 { 0, 0, 0 };
			if (status != initialisationStatus::initialisedAndApplied)
			{
				searchDensity(std::get<0>(particleSubstances_prt),
						particleSubstances_prt, concentrations, densities,
						normale,
						std::get<positionType::cell>(particlePosition_prt),
						rho1, r1, 1, boundVal, M);
				searchDensity(std::get<1>(particleSubstances_prt),
						particleSubstances_prt, concentrations, densities,
						normale,
						std::get<positionType::cell>(particlePosition_prt),
						rho2, r2, 1, boundVal, M);
			}

			for (std::size_t j = 0; j < particleVelocity().size(); ++j)
				velocity_arr[j] = particleVelocity()[j];
			for (std::size_t j = 0; j < particleDensityGrad().size(); ++j)
				gradRho_arr[j] = particleDensityGrad()[j];

			rho_arr[0] = rho1, rho_arr[1] = rho2;

			const auto cellIndex = std::get<positionType::cell>(
					particlePosition_prt);

			gradP_arr[0] = std::get<0>(gradP.cval()[cellIndex]()), gradP_arr[1] =
					std::get<1>(gradP.cval()[cellIndex]()), gradP_arr[2] =
					std::get<2>(gradP.cval()[cellIndex]());
			rho_arr[2] = densities[0].cval()[cellIndex];
			divU_arr[0] = divU.cval()[cellIndex];
			for (std::size_t i = 0; i < tensor::vsize; ++i)
				gradU_arr[i] = gradU.cval()[cellIndex]()[i];

			MPI_Send(velocity_arr.data(), vector::vsize, schemi_MPI_SCALAR,
					parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovVelocity),
					MPI_COMM_WORLD);
			MPI_Send(gradRho_arr.data(), vector::vsize, schemi_MPI_SCALAR,
					parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovGradRho),
					MPI_COMM_WORLD);
			MPI_Send(rho_arr.data(), 3, schemi_MPI_SCALAR, parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovRho),
					MPI_COMM_WORLD);
			MPI_Send(gradP_arr.data(), vector::vsize, schemi_MPI_SCALAR,
					parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovGradP),
					MPI_COMM_WORLD);
			MPI_Send(divU_arr.data(), 1, schemi_MPI_SCALAR, parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovDivU),
					MPI_COMM_WORLD);
			MPI_Send(gradU_arr.data(), tensor::vsize, schemi_MPI_SCALAR,
					parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovGradU),
					MPI_COMM_WORLD);
		}

		vector newPositionVector_prt;
		std::array<std::size_t, 1> nodeParticleLocatedNew { 0 };
		if (parallelism.isRoot())
		{
			MPI_Recv(velocity_arr.data(), vector::vsize, schemi_MPI_SCALAR,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovVelocity),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(gradRho_arr.data(), vector::vsize, schemi_MPI_SCALAR,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovGradRho),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(rho_arr.data(), 3, schemi_MPI_SCALAR,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovRho),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(gradP_arr.data(), vector::vsize, schemi_MPI_SCALAR,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovGradP),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(divU_arr.data(), 1, schemi_MPI_SCALAR,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovDivU),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(gradU_arr.data(), tensor::vsize, schemi_MPI_SCALAR,
					static_cast<int>(std::get<0>(nodeParticleLocatedOld)),
					static_cast<int>(MPIHandler::MPITags::GoncharovGradU),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			const auto particleVelocity = vector(velocity_arr[0],
					velocity_arr[1], velocity_arr[2]);
			const auto particleDensityGrad = vector(gradRho_arr[0],
					gradRho_arr[1], gradRho_arr[2]);

			const scalar rho1(rho_arr[0]), rho2(rho_arr[1]);

			const vector gradP_cell(gradP_arr[0], gradP_arr[1], gradP_arr[2]);
			const scalar rho_cell = rho_arr[2];
			const scalar divU_cell { divU_arr[0] };
			const tensor gradU_cell(gradU_arr[0], gradU_arr[1], gradU_arr[2],
					gradU_arr[3], gradU_arr[4], gradU_arr[5], gradU_arr[6],
					gradU_arr[7], gradU_arr[8]);

			particlesList[prt].timeIntegration(rho1, rho2, particleDensityGrad,
					particleVelocity, timestep, gradP_cell, rho_cell, divU_cell,
					gradU_cell);

			newPositionVector_prt = particlesList[prt].getPosition();
		}

		locateParticleNode(prt, newPositionVector_prt, nodeParticleLocatedNew);

		std::vector<std::size_t> newPositionData(0);
		std::vector<MPIHandler::mpi_scalar> newWeights_arr(0);
		if (parallelism.isRoot()
				|| (parallelism.mpi_rank == std::get<0>(nodeParticleLocatedNew)))
		{
			newPositionData.resize(2);
			newWeights_arr.resize(2);
		}

		if (parallelism.mpi_rank == std::get<0>(nodeParticleLocatedNew))
		{
			if (std::get<0>(nodeParticleLocatedNew)
					== std::get<0>(nodeParticleLocatedOld))
			{
				const auto indexes = findNewNearPosition(newPositionVector_prt,
						partLocOld_arr[0]);

				newPositionData[0] = std::get<0>(indexes);
				newPositionData[1] = std::get<1>(indexes);
			}
			else
			{
				newPositionData[0] = findNearCell(newPositionVector_prt);
				newPositionData[1] = findNearSurface(newPositionVector_prt,
						newPositionData[0]);
			}

			const std::array<std::size_t, 3> newParticlePos_prt {
					nodeParticleLocatedNew[0], newPositionData[0],
					newPositionData[1] };

			const auto newCellSurfWeight_prt = calculateWeights(
					newPositionVector_prt, newParticlePos_prt);

			for (std::size_t j = 0; j < newCellSurfWeight_prt.size(); ++j)
				newWeights_arr[j] = newCellSurfWeight_prt[j];

			MPI_Send(newPositionData.data(), 2, schemi_MPI_SIZE,
					parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
					MPI_COMM_WORLD);
			MPI_Send(newWeights_arr.data(), 2, schemi_MPI_SCALAR,
					parallelism.root,
					static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
					MPI_COMM_WORLD);
		}
		else if (parallelism.isRoot())
		{
			MPI_Recv(newPositionData.data(), 2, schemi_MPI_SIZE,
					static_cast<int>(std::get<0>(nodeParticleLocatedNew)),
					static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(newWeights_arr.data(), 2, schemi_MPI_SCALAR,
					static_cast<int>(std::get<0>(nodeParticleLocatedNew)),
					static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			std::get<positionType::cell>(particlePosition[prt]) =
					newPositionData[0];
			std::get<positionType::surface>(particlePosition[prt]) =
					newPositionData[1];

			cellSurfWeight[prt] = { newWeights_arr[0], newWeights_arr[1] };
		}
	}
#else
	for (std::size_t prt = 0; prt < particlesList.size(); ++prt)
	{
		const auto particleVelocity = calculateWeightedValue(uCell, uSurf,
				particlePosition[prt], cellSurfWeight[prt]);
		const auto particleDensityGrad = calculateWeightedValue(gradRhoCell,
				gradRhoSurf, particlePosition[prt], cellSurfWeight[prt]);

		const auto normale = particleDensityGrad / -particleDensityGrad.mag();

		scalar rho1 { 0 }, rho2 { 0 };
		vector r1 { 0, 0, 0 }, r2 { 0, 0, 0 };
		if (particleStatus[prt] != initialisationStatus::initialisedAndApplied)
		{
			searchDensity(std::get<0>(particlesList[prt].substances()),
					particlesList[prt].substances(), concentrations, densities,
					normale,
					std::get<positionType::cell>(particlePosition[prt]), rho1,
					r1, 1, boundVal, M);
			searchDensity(std::get<1>(particlesList[prt].substances()),
					particlesList[prt].substances(), concentrations, densities,
					normale,
					std::get<positionType::cell>(particlePosition[prt]), rho2,
					r2, 1, boundVal, M);
		}

		const auto cellIndex = std::get<positionType::cell>(
				particlePosition[prt]);

		particlesList[prt].timeIntegration(rho1, rho2, particleDensityGrad,
				particleVelocity, timestep, gradP.cval()[cellIndex],
				densities[0].cval()[cellIndex], divU.cval()[cellIndex],
				gradU.cval()[cellIndex]);

		const auto & newPositionVector(particlesList[prt].getPosition());

		if (!checkPosition(newPositionVector))
			throw exception("Particle out of domain.", errors::systemError);

		const auto indexes = findNewNearPosition(newPositionVector,
				std::get<positionType::cell>(particlePosition[prt]));

		std::get<positionType::cell>(particlePosition[prt]) = std::get<0>(
				indexes);
		std::get<positionType::surface>(particlePosition[prt]) = std::get<1>(
				indexes);

		cellSurfWeight[prt] = calculateWeights(newPositionVector,
				particlePosition[prt]);
	}
#endif
}

void schemi::instabilityParticlesHandler::writeOutput(
		const std::string & fieldDataDirectoryName, const scalar Time) const
{
	if (parallelism.isRoot() && modelUsed)
		for (std::size_t prt = 0; prt < particlesList.size(); ++prt)
		{
			std::string particleResultName { "./result/GoncharovModelParticle_" };
			particleResultName.append(std::to_string(prt));
			particleResultName.append(".tsv");

			std::string particleOutputName { fieldDataDirectoryName };
			particleOutputName.append("GoncharovModelParticleData_");
			particleOutputName.append(std::to_string(prt));
			particleOutputName.append(".dat");

			std::ofstream result { particleResultName, std::ios::app };

			if (!result.is_open())
				throw std::ofstream::failure(
						"File <<GoncharovModelParticle_>> is not opened.");

			result.precision(ioPrecision);

			const auto & prtcl = particlesList[prt];

			result << Time << '\t' << std::get<0>(prtcl.getPosition()()) << '\t'
					<< std::get<1>(prtcl.getPosition()()) << '\t'
					<< std::get<2>(prtcl.getPosition()()) << '\t'
					<< prtcl.getEta() << '\t'
					<< static_cast<int>(prtcl.getStatus()) << '\t'
					<< prtcl.getk0() << '\t' << prtcl.geteps0() << '\t'
					<< std::get<0>(prtcl.geta0()()) << '\t'
					<< std::get<1>(prtcl.geta0()()) << '\t'
					<< std::get<2>(prtcl.geta0()()) << '\t' << prtcl.getb0()
					<< '\t' << prtcl.getCurRe() << '\t'
					<< prtcl.getCurRelDeltaEta() << '\n';

			result.close();

			std::ofstream output { particleOutputName };

			if (!output.is_open())
				throw std::ofstream::failure(
						"File <<GoncharovModelParticleData_>> is not opened.");

			output.precision(ioPrecision);

			particlesList[prt].writeOutput(output);

			output.close();
		}
}

void schemi::instabilityParticlesHandler::checkTransitionToTurbulenceModel(
		const volumeField<scalar> & nuCell,
		const surfaceField<scalar> & nuSurface, volumeField<scalar> & k,
		volumeField<scalar> & epsilon, volumeField<vector> & a,
		volumeField<scalar> & b,
		const concentrationsPack<cubicCell> & concentrations,
		const boundaryConditionValue & boundVal, const scalar timestep)
{
	if (modelUsed)
	{
#ifdef MPI_VERSION
		for (std::size_t prt = 0; prt < listSize; ++prt)
		{
			std::vector<std::size_t> partLoc_arr(0);
			std::vector<MPIHandler::mpi_scalar> weights_arr(0), nodeData(0);

			std::array<std::size_t, 1> nodeParticleLocated { 0 };

			if (parallelism.isRoot())
				std::get<0>(nodeParticleLocated) = std::get<positionType::rank>(
						particlePosition[prt]);

			MPI_Bcast(nodeParticleLocated.data(), 1, schemi_MPI_SIZE,
					parallelism.root,
					MPI_COMM_WORLD);

			if (parallelism.isRoot()
					|| (parallelism.mpi_rank == std::get<0>(nodeParticleLocated)))
			{
				partLoc_arr.resize(2);
				weights_arr.resize(2);
				nodeData.resize(7);
			}

			if (parallelism.isRoot())
			{
				for (std::size_t j = 0; j < 2; ++j)
				{
					partLoc_arr[j] = particlePosition[prt][j + 1];
					weights_arr[j] = cellSurfWeight[prt][j];
				}

				MPI_Send(partLoc_arr.data(), 2, schemi_MPI_SIZE,
						static_cast<int>(std::get<0>(nodeParticleLocated)),
						static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
						MPI_COMM_WORLD);
				MPI_Send(weights_arr.data(), 2, schemi_MPI_SCALAR,
						static_cast<int>(std::get<0>(nodeParticleLocated)),
						static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
						MPI_COMM_WORLD);
			}
			else if (parallelism.mpi_rank == std::get<0>(nodeParticleLocated))
			{
				MPI_Recv(partLoc_arr.data(), 2, schemi_MPI_SIZE,
						parallelism.root,
						static_cast<int>(MPIHandler::MPITags::GoncharovPosition),
						MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				MPI_Recv(weights_arr.data(), 2, schemi_MPI_SCALAR,
						parallelism.root,
						static_cast<int>(MPIHandler::MPITags::GoncharovWeights),
						MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

				const std::array<std::size_t, 3> particlePosition_prt {
						nodeParticleLocated[0], partLoc_arr[0], partLoc_arr[1] };
				const std::array<scalar, 2> cellSurfWeight_prt { weights_arr[0],
						weights_arr[1] };

				const auto nu = calculateWeightedValue(nuCell, nuSurface,
						particlePosition_prt, cellSurfWeight_prt);

				const auto cellRadius(
						meshRef.cells()[std::get<positionType::cell>(
								particlePosition_prt)].rC());
				const auto surfaceRadius(
						meshRef.surfaces()[std::get<positionType::surface>(
								particlePosition_prt)].rC());

				nodeData[0] = nu;
				for (std::size_t i = 0; i < vector::vsize; ++i)
				{
					nodeData[i + 1] = cellRadius()[i];
					nodeData[i + 4] = surfaceRadius()[i];
				}

				MPI_Send(nodeData.data(), 7, schemi_MPI_SCALAR,
						parallelism.root,
						static_cast<int>(MPIHandler::MPITags::GoncharovViscosityAndRadiuses),
						MPI_COMM_WORLD);
			}

			if (parallelism.isRoot())
			{
				MPI_Recv(nodeData.data(), 7, schemi_MPI_SCALAR,
						static_cast<int>(std::get<0>(nodeParticleLocated)),
						static_cast<int>(MPIHandler::MPITags::GoncharovViscosityAndRadiuses),
						MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

				const scalar nu = nodeData[0];
				const vector cellRadius = { nodeData[1], nodeData[2],
						nodeData[3] };
				const vector surfaceRadius = { nodeData[4], nodeData[5],
						nodeData[6] };

				if (particleStatus[prt] == initialisationStatus::notInitialised)
				{
					const auto particleTransitionStatus =
							particlesList[prt].checkTransition(nu, timestep,
									cellRadius, surfaceRadius);

					if (particleTransitionStatus
							== interfaceStatus::developedResolvable)
						particleStatus[prt] =
								initialisationStatus::initialisedButNotApplied;
				}
			}
		}

		for (std::size_t prt = 0; prt < listSize; ++prt)
		{
			std::array<int, 1> partStatusCast;

			if (parallelism.isRoot())
				std::get<0>(partStatusCast) =
						static_cast<int>(particleStatus[prt]);

			MPI_Bcast(partStatusCast.data(), 1, MPI_INT, parallelism.root,
			MPI_COMM_WORLD);

			const initialisationStatus partStatus { std::get<0>(partStatusCast) };

			if (partStatus == initialisationStatus::initialisedButNotApplied)
			{
				std::array<std::size_t, 2> substIndexes;
				std::array<MPIHandler::mpi_scalar, 4> interfPos;
				std::array<MPIHandler::mpi_scalar, 6> turbInitData;

				if (parallelism.isRoot())
				{
					const auto & particle_prt = particlesList[prt];

					substIndexes = particle_prt.substances();

					interfPos = { std::get<0>(particle_prt.getPosition()()),
							std::get<1>(particle_prt.getPosition()()), std::get<
									2>(particle_prt.getPosition()()),
							particle_prt.getRadiusOfInfluence() };

					turbInitData = { particle_prt.getk0(),
							particle_prt.geteps0(), std::get<0>(
									particle_prt.geta0()()), std::get<1>(
									particle_prt.geta0()()), std::get<2>(
									particle_prt.geta0()()),
							particle_prt.getb0() };
				}

				MPI_Bcast(substIndexes.data(), 2, schemi_MPI_SIZE,
						parallelism.root,
						MPI_COMM_WORLD);
				MPI_Bcast(interfPos.data(), 4, schemi_MPI_SCALAR,
						parallelism.root,
						MPI_COMM_WORLD);
				MPI_Bcast(turbInitData.data(), 6, schemi_MPI_SCALAR,
						parallelism.root,
						MPI_COMM_WORLD);

				const auto profile = formProfile(substIndexes,
						vector(std::get<0>(interfPos), std::get<1>(interfPos),
								std::get<2>(interfPos)), std::get<3>(interfPos),
						concentrations, boundVal);

				const auto k0 = std::get<0>(turbInitData);
				const auto eps0 = std::get<1>(turbInitData);
				const vector a0 = { std::get<2>(turbInitData), std::get<3>(
						turbInitData), std::get<4>(turbInitData) };
				const auto b0 = std::get<5>(turbInitData);

				for (std::size_t i = 0; i < profile.size(); ++i)
				{
					const scalar profileFactor = profile.cval()[i];

					k.val()[i] += profileFactor * k0;
					epsilon.val()[i] += profileFactor * eps0;
					a.val()[i] += profileFactor * a0;
					b.val()[i] += profileFactor * b0;
				}

				if (parallelism.isRoot())
					particleStatus[prt] =
							initialisationStatus::initialisedAndApplied;
			}
		}
#else
		for (std::size_t prt = 0; prt < particlesList.size(); ++prt)
		{
			const auto nu = calculateWeightedValue(nuCell, nuSurface,
					particlePosition[prt], cellSurfWeight[prt]);

			const auto & cellRadius(
					meshRef.cells()[std::get<positionType::cell>(
							particlePosition[prt])].rC());
			const auto & surfaceRadius(
					meshRef.surfaces()[std::get<positionType::surface>(
							particlePosition[prt])].rC());

			if (particleStatus[prt] == initialisationStatus::notInitialised)
			{
				const auto particleTransitionStatus =
						particlesList[prt].checkTransition(nu, timestep,
								cellRadius, surfaceRadius);

				if (particleTransitionStatus
						== interfaceStatus::developedResolvable)
					particleStatus[prt] =
							initialisationStatus::initialisedButNotApplied;
			}
		}

		for (std::size_t prt = 0; prt < particlesList.size(); ++prt)
			if (particleStatus[prt]
					== initialisationStatus::initialisedButNotApplied)
			{
				const auto & particle_prt = particlesList[prt];

				const auto profile = formProfile(particle_prt.substances(),
						particle_prt.getPosition(),
						particle_prt.getRadiusOfInfluence(), concentrations,
						boundVal);

				for (std::size_t i = 0; i < profile.size(); ++i)
				{
					const scalar profileFactor = profile.cval()[i];

					k.val()[i] += profileFactor * particle_prt.getk0();
					epsilon.val()[i] += profileFactor * particle_prt.geteps0();
					a.val()[i] += profileFactor * particle_prt.geta0();
					b.val()[i] += profileFactor * particle_prt.getb0();
				}

				particleStatus[prt] =
						initialisationStatus::initialisedAndApplied;
			}
#endif
	}
}

schemi::volumeField<schemi::scalar> schemi::instabilityParticlesHandler::generationArea(
		const concentrationsPack<cubicCell> & concentrations) const noexcept
{
	volumeField<scalar> generationFlag(concentrations.v[0].meshRef(), 0.0);

	if (!modelUsed)
	{
		generationFlag.val() = 1;
		return generationFlag;
	}

	constexpr scalar purityParameterExpanded = 1 - (1 - purityParameter) * 0.1;

#ifdef MPI_VERSION
	for (std::size_t prt = 0; prt < listSize; ++prt)
	{
		std::array<int, 1> partStatusCast;
		std::array<std::size_t, 2> substancesCast;
		std::array<MPIHandler::mpi_scalar, 4> positionCast;

		if (parallelism.isRoot())
		{
			std::get<0>(partStatusCast) = static_cast<int>(particleStatus[prt]);
			substancesCast = particlesList[prt].substances();
			const auto & interfPos = particlesList[prt].getPosition();
			positionCast = { std::get<0>(interfPos()), std::get<1>(interfPos()),
					std::get<2>(interfPos()),
					particlesList[prt].getRadiusOfInfluence() };
		}

		MPI_Bcast(partStatusCast.data(), 1, MPI_INT, parallelism.root,
		MPI_COMM_WORLD);

		const initialisationStatus partStatus { std::get<0>(partStatusCast) };

		if (partStatus == initialisationStatus::initialisedAndApplied)
		{
			MPI_Bcast(substancesCast.data(), 2, schemi_MPI_SIZE,
					parallelism.root,
					MPI_COMM_WORLD);
			MPI_Bcast(positionCast.data(), 4, schemi_MPI_SCALAR,
					parallelism.root,
					MPI_COMM_WORLD);

			const auto & mixSubs = substancesCast;
			const auto interfacePosition = vector(std::get<0>(positionCast),
					std::get<1>(positionCast), std::get<2>(positionCast));
			const scalar radiusOfInfluence = std::get<3>(positionCast);

			for (std::size_t i = 0; i < concentrations.v[0].size(); ++i)
			{
				const auto x1 =
						concentrations.v[std::get<0>(mixSubs) + 1].cval()[i]
								/ concentrations.v[0].cval()[i];
				const auto x2 =
						concentrations.v[std::get<1>(mixSubs) + 1].cval()[i]
								/ concentrations.v[0].cval()[i];

				const auto deltaR = meshRef.cells()[i].rC() - interfacePosition;

				if ((x1 <= purityParameterExpanded)
						&& (x2 <= purityParameterExpanded)
						&& (deltaR.mag() < radiusOfInfluence))
					generationFlag.val()[i] = 1;
			}
		}
	}
#else
	for (std::size_t prt = 0; prt < particlesList.size(); ++prt)
		if (particleStatus[prt] == initialisationStatus::initialisedAndApplied)
		{
			const auto & mixSubs = particlesList[prt].substances();
			const auto & interfacePosition = particlesList[prt].getPosition();

			for (std::size_t i = 0; i < concentrations.v[0].size(); ++i)
			{
				const auto x1 =
						concentrations.v[std::get<0>(mixSubs) + 1].cval()[i]
								/ concentrations.v[0].cval()[i];
				const auto x2 =
						concentrations.v[std::get<1>(mixSubs) + 1].cval()[i]
								/ concentrations.v[0].cval()[i];

				const auto deltaR = meshRef.cells()[i].rC() - interfacePosition;

				if ((x1 <= purityParameterExpanded)
						&& (x2 <= purityParameterExpanded)
						&& (deltaR.mag()
								< particlesList[prt].getRadiusOfInfluence()))
					generationFlag.val()[i] = 1;
			}
		}
#endif

	return generationFlag;
}
