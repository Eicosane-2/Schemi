/*
 * main.cpp
 *
 *  Created on: 2019/11/12
 *      Author: Maxim Boldyrev
 */

#include <cstdlib>
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "timestepEnum.hpp"
#include "chemicalReactionsEnum.hpp"

#include "abstractChemicalKinetics.hpp"
#include "abstractFlowSolver.hpp"
#include "abstractLimiter.hpp"
#include "abstractStepSolver.hpp"
#include "fabricFunctions.hpp"
#include "inputString.hpp"
#include "MPIHandler.hpp"
#include "output.hpp"
#include "phaseInitialization.hpp"
#include "scalar.hpp"
#include "secondOrderStepSolver.hpp"
#include "structForOutput.hpp"
#include "thirdOrderStepSolver.hpp"

#ifdef MPI_VERSION
int main(int argc, char * argv[])
#else
int main()
#endif
{
	using namespace schemi;

	try
	{

		const auto startTime = std::chrono::high_resolution_clock::now();
		scalar timeForTVD { 0 }, timeForHancock { 0 }, timeForFlowCalculation {
				0 }, timeForTimeIntegration { 0 }, timeForDiffusion { 0 };

		int rank(0), size(1);
#ifdef MPI_VERSION
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

		MPIHandler parallelism(rank, size);

#ifdef MPI_VERSION
		if (parallelism.mpi_size <= 1)
			throw exception(
					"Program's work with number of nodes less than 2 has not been tested.",
					errorsEnum::initializationError);
#endif

		vector systemSize;
		scalar timeOfCalculation, minTime, Courant, Time { 0 };
		std::pair<scalar, scalar> timestepCoeffs { 0.25, 0.5 };
		std::size_t numberOfCells_x, numberOfCells_y, numberOfCells_z,
				frequencyOfOutput, frequencyOfOutputWidth, numberOfComponents,
				numberOfZones, numberOfIterations, nsteps { 0 }, nouts { 0 },
				noutsW { 0 };
		std::string universalGasConstant, equationOfState,
				typeOfTVDLimeterString, diffusionONString, turbulenceONString,
				matrixSolverString, gravitationONString, sourceTypeString,
				flowSolwerString, dimensionsOfTask, sourceTimeFlagString,
				linearFlagString, thirdOrderString, readFromOutput;
		bool diffusionFlag, gravitationFlag, linearFlag;
		vector g { 0 }, gDelta { 0 };
		dimensionsEnum dimensionsFlag;
		timestepEnum sourceTimeFlag;
		std::pair<bool, scalar> constTimeStep { false, 0. };

		std::string skipBuffer;

		/*Reading gravitational acceleration.*/
		{
			std::ifstream gFile { "./set/g.txt" };
			if (gFile.is_open())
				std::cout << "./set/g.txt is opened." << std::endl;
			else
				throw exception("./set/g.txt not found.",
						errorsEnum::initializationError);

			gFile >> skipBuffer >> gravitationONString >> skipBuffer
					>> g.v_r()[0] >> g.v_r()[1] >> g.v_r()[2] >> skipBuffer
					>> gDelta.v_r()[0] >> gDelta.v_r()[1] >> gDelta.v_r()[2];
		}

		/*Reading initial and boundary conditions.*/
		{
			std::ifstream mainParametersFile { "./set/main.txt" };
			if (mainParametersFile.is_open())
				std::cout << "./set/main.txt is opened." << std::endl;
			else
				throw exception("./set/main.txt not found.",
						errorsEnum::initializationError);

			mainParametersFile

			>> skipBuffer >> systemSize.v_r()[0] >> systemSize.v_r()[1]
					>> systemSize.v_r()[2]

					>> skipBuffer >> numberOfCells_x >> numberOfCells_y
					>> numberOfCells_z

					>> skipBuffer >> timeOfCalculation

					>> skipBuffer >> Courant

					>> skipBuffer >> frequencyOfOutput

					>> skipBuffer >> minTime

					>> skipBuffer >> numberOfComponents

					>> skipBuffer >> numberOfZones

					>> skipBuffer >> universalGasConstant

					>> skipBuffer >> equationOfState

					>> skipBuffer >> typeOfTVDLimeterString

					>> skipBuffer >> numberOfIterations

					>> skipBuffer >> diffusionONString

					>> skipBuffer >> turbulenceONString

					>> skipBuffer >> matrixSolverString

					>> skipBuffer >> sourceTypeString

					>> skipBuffer >> flowSolwerString

					>> skipBuffer >> dimensionsOfTask

					>> skipBuffer >> constTimeStep.second

					>> skipBuffer >> sourceTimeFlagString

					>> skipBuffer >> linearFlagString

					>> skipBuffer >> thirdOrderString

					>> skipBuffer >> readFromOutput;

			mainParametersFile.close();
		}

		/*Setting flags.*/
		if (linearFlagString == "linear")
			linearFlag = true;
		else if (linearFlagString == "star")
			linearFlag = false;
		else
			throw exception("Unknown reconstruction flag.",
					errorsEnum::initializationError);

		frequencyOfOutputWidth = std::min(frequencyOfOutput, std::size_t(1000));

		if (constTimeStep.second > 0.)
			constTimeStep.first = true;

		{
			const auto eqCount = std::count(sourceTimeFlagString.begin(),
					sourceTimeFlagString.end(), '=');
			const auto comCount = std::count(sourceTimeFlagString.begin(),
					sourceTimeFlagString.end(), ',');

			const auto eqFind = sourceTimeFlagString.find('=');
			const auto comFind = sourceTimeFlagString.find(',');

			if (((eqCount != 1) && (comCount != 1)) || (eqFind > comFind))
				throw exception("Wrong format of time-step parameters settings",
						errorsEnum::initializationError);

			std::string sourceTimeFlagType;
			auto eqIt = sourceTimeFlagString.begin();
			auto comIt = sourceTimeFlagString.begin();
			for (auto STF_striter = sourceTimeFlagString.begin();
					STF_striter != sourceTimeFlagString.end(); ++STF_striter)
			{
				if (*STF_striter == '=')
					eqIt = STF_striter;
				else if (*STF_striter == ',')
				{
					comIt = STF_striter;

					break;
				}
			}

			{
				sourceTimeFlagType = std::string(sourceTimeFlagString.begin(),
						eqIt);

				auto num1_beg = eqIt;
				num1_beg++;

				auto num2_beg = comIt;
				num2_beg++;

				const std::string number1(num1_beg, comIt);
				const std::string number2(num2_beg, sourceTimeFlagString.end());

				timestepCoeffs.first = std::stod(number1);
				timestepCoeffs.second = std::stod(number2);

				std::cout << "Timestep type: " << sourceTimeFlagType
						<< std::endl;
				std::cout << "Source and diffusion timestep coefficients: "
						<< timestepCoeffs.first << ' ' << timestepCoeffs.second
						<< std::endl;
			}

			if (sourceTimeFlagType == "Courant")
			{
				std::cout << "Courant time-step on." << std::endl;

				sourceTimeFlag = timestepEnum::CourantTimeStep;
			}
			else if (sourceTimeFlagType == "Courant-Source")
			{
				std::cout << "Courant and source time-step on." << std::endl;

				sourceTimeFlag = timestepEnum::CourantAndSourceTimeStep;
			}
			else if (sourceTimeFlagType == "Courant-Source-Diffusion")
			{
				std::cout << "Courant, source and diffusion time-step on."
						<< std::endl;

				sourceTimeFlag =
						timestepEnum::CourantAndSourceAndDiffusionTimeStep;
			}
			else
				throw exception("Unknown source time-step flag.",
						errorsEnum::initializationError);
		}

		if (gravitationONString == "on")
			gravitationFlag = true;
		else if (gravitationONString == "off")
			gravitationFlag = false;
		else
			throw exception("Unknown gravitation flag.",
					errorsEnum::initializationError);

		if (diffusionONString == "on")
			diffusionFlag = true;
		else if (diffusionONString == "off")
			diffusionFlag = false;
		else
			throw exception("Unknown diffusion flag.",
					errorsEnum::initializationError);

		std::unique_ptr<abstractLimiter> limiter(
				createLimiter(typeOfTVDLimeterString));

		auto [msolver, msolverEnthFl] = createMatrixSolver(matrixSolverString,
				dimensionsOfTask, numberOfIterations);

		std::unique_ptr<abstractFlowSolver> fsolver(
				createFlowSolver(flowSolwerString));

		if (dimensionsOfTask == "1D")
			dimensionsFlag = dimensionsEnum::task1D;
		else if (dimensionsOfTask == "2D")
			dimensionsFlag = dimensionsEnum::task2D;
		else if (dimensionsOfTask == "3D")
			dimensionsFlag = dimensionsEnum::task3D;
		else
			throw exception("Unknown dimensions flag.",
					errorsEnum::initializationError);

		bool thirdOrder;
		if (thirdOrderString == "SecondOrder")
			thirdOrder = false;
		else if (thirdOrderString == "ThirdOrder")
			thirdOrder = true;
		else
			throw exception("Unknown gas dynamics approximation order.",
					errorsEnum::initializationError);

		std::size_t readDataPoint;
		if ((readFromOutput == "no") || (readFromOutput == "0"))
		{
			readDataPoint = 0;

			if (parallelism.isRoot())
			{
				const std::string resultFileName { "./result" };

				if (!std::filesystem::exists(resultFileName))
					std::filesystem::create_directory(resultFileName);

				std::cout << "Time = " << Time << std::endl;
				{
					/*Recreation of Time.tsv, if it exist.*/
					std::string timeFileName("./result/Time.tsv");
					std::ofstream timeFile(timeFileName);
					timeFile.close();
				}
				if (dimensionsFlag == dimensionsEnum::task1D)
				{
					/*Recreation of timeWidth.tsv, if it exist.*/
					std::string timeWidthFileName("./result/timeWidth.tsv");
					std::ofstream timeWidthFile(timeWidthFileName);
					timeWidthFile << "Time" << '\t' << "Width" << std::endl;
					timeWidthFile.close();
				}
			}
		}
		else
		{
			readDataPoint = std::stoul(readFromOutput);

			if (parallelism.isRoot())
			{
				const std::string resultFileName { "./result" };
				if (!std::filesystem::exists(resultFileName))
					throw exception(
							std::string("<<result>> directory doesn't exist."),
							errorsEnum::systemError);

				std::string timeFileName("./result/Time.tsv");
				std::ifstream timeFile(timeFileName);
				if (!timeFile.is_open())
					throw exception(std::string("Couldn't open Time.tsv"),
							errorsEnum::systemError);
				timeFile.precision(ioPrecision);

				std::string timeString;
				std::size_t lineNumber { 0 };
				std::size_t lineNumberEnd { 0 };
				bool isLastOutput { false };
				while (std::getline(timeFile, skipBuffer))
				{
					lineNumberEnd++;

					if ((lineNumberEnd - 1) == readDataPoint)
					{
						lineNumber = lineNumberEnd;

						timeString = skipBuffer;
					}
				}

				if (lineNumber == 0)
					throw exception(
							std::string(
									"Appropriate line wasn't found in Time.tsv"),
							errorsEnum::systemError);

				if (lineNumber == 1)
					throw exception(
							std::string("Time.tsv contains only one line."),
							errorsEnum::systemError);

				if (lineNumberEnd == lineNumber)
					isLastOutput = true;

				auto tab = timeString.begin();

				for (auto ts_striter = timeString.begin();
						ts_striter != timeString.end(); ++ts_striter)
					if (*ts_striter == '\t')
					{
						tab = ts_striter;
						break;
					}

				std::string noutsString = std::string(timeString.begin(), tab);

				auto TimeIter = tab;
				TimeIter++;

				std::string TimeString = std::string(TimeIter,
						timeString.end());

				nouts = std::stoul(noutsString);
				Time = std::stod(TimeString);

				const auto remainingTime { timeOfCalculation - Time };

				if (remainingTime <= 0)
					throw exception(
							"Time of calculation lesser or equal than time of executed calculation.",
							errorsEnum::initializationError);

				if (nouts != readDataPoint)
					throw exception(
							std::string(
									"Number of output does not concur with <<readDataPoint>>."),
							errorsEnum::systemError);

				std::vector<std::string> savedStrings(lineNumber);
				timeFile.clear();
				timeFile.seekg(0, timeFile.beg);
				for (std::size_t str_i = 0; str_i < lineNumber; ++str_i)
					std::getline(timeFile, savedStrings[str_i]);

				timeFile.close();

				/*Recreation of Time.tsv.*/
				std::ofstream timeFileNew(timeFileName);
				if (!timeFileNew.is_open())
					throw exception(std::string("Couldn't open Time.tsv"),
							errorsEnum::systemError);
				timeFileNew.precision(ioPrecision);

				for (std::size_t str_i = 0; str_i < lineNumber; ++str_i)
					timeFileNew << savedStrings[str_i] << '\n';

				timeFileNew.close();

				if (isLastOutput)
				{
					std::ifstream numberOfStepsFile("./timeOfCalculation.tsv");
					if (!numberOfStepsFile.is_open())
						throw exception(
								std::string(
										"Couldn't open timeOfCalculation.tsv"),
								errorsEnum::systemError);

					while (!numberOfStepsFile.eof())
					{
						std::string word;
						numberOfStepsFile >> word;

						if (word == "steps")
						{
							numberOfStepsFile >> word;

							nsteps = std::stoul(word);

							break;
						}
					}
				}
				else
					nsteps = nouts * frequencyOfOutput;

				if (dimensionsFlag == dimensionsEnum::task1D)
				{
					std::string timeWidthFileName("./result/timeWidth.tsv");
					std::ifstream timeWidthFile(timeWidthFileName);
					if (!timeWidthFile.is_open())
						throw exception(
								std::string("Couldn't open timeWidth.tsv"),
								errorsEnum::systemError);

					lineNumber = 0;
					std::string widthData;
					while (std::getline(timeWidthFile, widthData))
					{
						lineNumber++;

						if (lineNumber == 1)
							continue;
						else
						{
							scalar timeLastWidth;

							tab = widthData.begin();

							for (auto ts_striter = widthData.begin();
									ts_striter != widthData.end(); ++ts_striter)
								if (*ts_striter == '\t')
								{
									tab = ts_striter;
									break;
								}

							timeLastWidth = std::stod(
									std::string(widthData.begin(), tab));

							if (timeLastWidth == Time)
								break;
							else if (timeLastWidth > Time)
							{
								lineNumber--;
								break;
							}
							else
								continue;
						}
					}

					noutsW = lineNumber - 1;

					savedStrings.resize(lineNumber);
					timeWidthFile.clear();
					timeWidthFile.seekg(0, timeFile.beg);
					for (std::size_t str_i = 0; str_i < lineNumber; ++str_i)
						std::getline(timeWidthFile, savedStrings[str_i]);

					timeWidthFile.close();

					/*Recreation of timeWidth.tsv.*/
					std::ofstream timeWidthFileNew(timeWidthFileName);
					if (!timeWidthFileNew.is_open())
						throw exception(
								std::string("Couldn't open timeWidth.tsv"),
								errorsEnum::systemError);
					timeWidthFileNew.precision(ioPrecision);

					for (std::size_t str_i = 0; str_i < lineNumber; ++str_i)
						timeWidthFileNew << savedStrings[str_i] << '\n';

					timeWidthFileNew.close();
				}
			}

#ifdef MPI_VERSION
			{
				std::size_t noutsBroadcast[1], nstepsBroadcast[1],
						noutsWBroadcast[1];
				MPIHandler::mpi_scalar TimeBroadcast[1];

				if (parallelism.isRoot())
				{
					noutsBroadcast[0] = nouts;
					nstepsBroadcast[0] = nsteps;
					noutsWBroadcast[0] = noutsW;
					TimeBroadcast[0] = Time;
				}

				MPI_Bcast(noutsBroadcast, 1, schemi_MPI_SIZE, parallelism.root,
				MPI_COMM_WORLD);
				MPI_Bcast(nstepsBroadcast, 1, schemi_MPI_SIZE, parallelism.root,
				MPI_COMM_WORLD);
				MPI_Bcast(noutsWBroadcast, 1, schemi_MPI_SIZE, parallelism.root,
				MPI_COMM_WORLD);
				MPI_Bcast(TimeBroadcast, 1, MPI_DOUBLE, parallelism.root,
				MPI_COMM_WORLD);

				nouts = noutsBroadcast[0];
				nsteps = nstepsBroadcast[0];
				noutsW = noutsWBroadcast[0];
				Time = TimeBroadcast[0];
			}
#endif
		}

		/*Read common boundary condition.*/
		std::vector<boundaryConditionType> commonConditions { 6,
				boundaryConditionType::blank };
		{
			switch (dimensionsFlag)
			{
			case dimensionsEnum::task3D:
				commonConditions[2] = boundaryConditionType::calculated;
				commonConditions[5] = boundaryConditionType::calculated;
				[[fallthrough]];
			case dimensionsEnum::task2D:
				commonConditions[3] = boundaryConditionType::calculated;
				commonConditions[4] = boundaryConditionType::calculated;
				[[fallthrough]];
			case dimensionsEnum::task1D:
			default:
				commonConditions[0] = boundaryConditionType::calculated;
				commonConditions[1] = boundaryConditionType::calculated;
				break;
			}
		}

		/*Create mesh.*/
		mesh * const meshObjPointer = mesh::instance(); /*Self-destructing singleton.*/
		{
			auto parallelNodeSystemSize =
					parallelism.correctParallelepipedVector(systemSize);

			switch (dimensionsFlag)
			{
			case dimensionsEnum::task3D:
				meshObjPointer->threeDParallelepiped(parallelNodeSystemSize,
						numberOfCells_x, numberOfCells_y, numberOfCells_z,
						commonConditions);
				break;
			case dimensionsEnum::task2D:
				meshObjPointer->twoDParallelepiped(parallelNodeSystemSize,
						numberOfCells_x, numberOfCells_y, commonConditions);
				break;
			case dimensionsEnum::task1D:
			default:
				meshObjPointer->oneDParallelepiped(parallelNodeSystemSize,
						numberOfCells_x, commonConditions);
				break;
			}
		}

		auto & mesh = *meshObjPointer;

		parallelism.initialiseBuffersSize(mesh);
		parallelism.initializeParalleMeshData(mesh);

		/*Check number of components.*/
		if (numberOfComponents > 9)
			throw exception("More than 9 components.",
					errorsEnum::initializationError);

		/*Creating fields.*/
		auto [gasPhase, enthalpyFlowFlag, molMassDiffusionFlag] =
				phaseInitialization(numberOfComponents, numberOfZones, mesh,
						commonConditions, parallelism, turbulenceONString,
						sourceTypeString, universalGasConstant, equationOfState,
						readDataPoint);

		volumeField<scalar> sonicSpeed { mesh, 0 };
		sonicSpeed.ref_r() = std::sqrt(
				gasPhase->phaseThermodynamics->sqSonicSpeed(
						gasPhase->concentration.p, gasPhase->density[0].ref(),
						gasPhase->internalEnergy.ref(),
						gasPhase->pressure.ref()));

		/*Calculate effective length for time-step calculation and set first time-step.*/
		volumeField<scalar> minimalLengthScale { mesh, 0 };
		for (std::size_t i = 0; i < mesh.cellsSize(); ++i)
		{
			const auto edge1 =
					(mesh.cells()[i].rX00() - mesh.cells()[i].r000()).mag();
			const auto edge2 =
					(mesh.cells()[i].r0Y0() - mesh.cells()[i].r000()).mag();
			const auto edge3 =
					(mesh.cells()[i].r00Z() - mesh.cells()[i].r000()).mag();

			const auto minEdge = std::min(std::min(edge1, edge2), edge3);

			minimalLengthScale.ref_r()[i] = 1. / minEdge;
		}

		/*Set small initial step for safe source integration.*/
		if (!constTimeStep.first)
		{
			mesh.setTimestep(minTime);
			mesh.setTimestep(std::min(mesh.timestep(), timeOfCalculation));
		}
		else
			mesh.setTimestep(constTimeStep.second);

		chemicalReactionsEnum chemReactionFlag;
		{
			std::ifstream chem { "./set/chemicalKinetics.txt" };

			if (chem.is_open())
				std::cout << "./set/chemicalKinetics.txt is opened."
						<< std::endl;
			else
				throw exception("./set/chemicalKinetics.txt not found.",
						errorsEnum::initializationError);

			std::string reactionName;

			chem >> skipBuffer >> reactionName;

			if (reactionName == "no")
				chemReactionFlag = chemicalReactionsEnum::noReaction;
			else if (reactionName == "Cl2")
				chemReactionFlag = chemicalReactionsEnum::Cl2Dissociation;
			else if (reactionName == "Cl2H2")
				chemReactionFlag = chemicalReactionsEnum::Cl2H2Dissociation;
			else if (reactionName == "H2Cl2Combustion")
				chemReactionFlag = chemicalReactionsEnum::H2Cl2Combustion;
			else if (reactionName == "NO2Disproportionation")
				chemReactionFlag = chemicalReactionsEnum::NO2Disproportionation;
			else
				throw exception("Unknown chemical reaction model.",
						errorsEnum::initializationError);
		}

		std::unique_ptr<abstractChemicalKinetics> chmk(
				createChemicalKinetics(*gasPhase, chemReactionFlag));

		/*Write initial conditions.*/
		{
			if (!readDataPoint)
			{
				structForOutput outputData(parallelism, mesh,
						numberOfComponents);
				if (parallelism.isRoot())
					outputData.setSizes();

				outputData.collectParallelData(*gasPhase, gasPhase->tNu,
						sonicSpeed);
				if (parallelism.isRoot())
				{
					output::dataOutput(outputData, nouts, Time);

					if (dimensionsFlag == dimensionsEnum::task1D)
						output::mixedZoneWidth1D(outputData, Time);
				}
			}
			nouts++, noutsW++;
		}

		boundaryConditionValue boundaryConditionValueCalc(
				*gasPhase->turbulenceSources->turbPar, *gasPhase,
				*(gasPhase->phaseThermodynamics));

		skipBuffer.clear();

		parallelism.correctBoundaryValues(*gasPhase);

		std::unique_ptr<abstractStepSolver> stepSolver;
		if (thirdOrder)
			stepSolver = std::make_unique<thirdOrderStepSolver>(*gasPhase,
					*limiter, *fsolver, gravitationFlag, g,
					boundaryConditionValueCalc, timeForTVD, timeForHancock,
					timeForFlowCalculation, timeForTimeIntegration, parallelism,
					diffusionFlag, *msolver, *msolverEnthFl, timestepCoeffs,
					timeForDiffusion, commonConditions, enthalpyFlowFlag,
					linearFlag, boundaryConditionValueCalc, minimalLengthScale,
					sourceTimeFlag, molMassDiffusionFlag, *chmk);
		else
			stepSolver = std::make_unique<secondOrderStepSolver>(*gasPhase,
					*limiter, *fsolver, gravitationFlag, g,
					boundaryConditionValueCalc, timeForTVD, timeForHancock,
					timeForFlowCalculation, timeForTimeIntegration, parallelism,
					diffusionFlag, *msolver, *msolverEnthFl, timestepCoeffs,
					timeForDiffusion, commonConditions, enthalpyFlowFlag,
					linearFlag, boundaryConditionValueCalc, minimalLengthScale,
					sourceTimeFlag, molMassDiffusionFlag, *chmk);

#ifdef MPI_VERSION
		(MPI_Barrier(MPI_COMM_WORLD));
#endif

		do
		{
			nsteps++;

			try
			{
				stepSolver->calculateStep();
			} catch (const std::exception & exception)
			{
				std::clog << exception.what() << std::endl;
				std::clog << "Step number " << nsteps << '.' << std::endl;
				std::clog << "Time-step " << mesh.timestep() << '.'
						<< std::endl;
#ifndef MPI_VERSION
				inputString();
				std::clog << "Emergency output." << std::endl;
				{
					structForOutput outputData(parallelism, mesh,
							numberOfComponents);
					if (parallelism.isRoot())
						outputData.setSizes();

					outputData.collectParallelData(*gasPhase, gasPhase->tNu,
							sonicSpeed);
					if (parallelism.isRoot())
						output::dataOutput(outputData, nouts, Time);

					nouts++;
				}
				inputString();
#endif
				std::exit(EXIT_FAILURE);
			}

			if (gravitationFlag)
				g += gDelta * mesh.timestep();

			Time += mesh.timestep();

			/*Calculation of new time-step through speed of sound.*/
			if (!constTimeStep.first)
			{
				sonicSpeed.ref_r() = std::sqrt(
						gasPhase->phaseThermodynamics->sqSonicSpeed(
								gasPhase->concentration.p,
								gasPhase->density[0].ref(),
								gasPhase->internalEnergy.ref(),
								gasPhase->pressure.ref()));

				std::valarray<scalar> signalSpeed(
						gasPhase->velocity.ref().size());

				for (std::size_t i = 0; i < mesh.cellsSize(); ++i)
					signalSpeed[i] = gasPhase->velocity.ref()[i].mag()
							+ sonicSpeed.ref()[i];

				signalSpeed *= minimalLengthScale.ref();

				const scalar maxSignalSpeed { signalSpeed.max() };

				mesh.setTimestep(Courant / maxSignalSpeed);
				if (sourceTimeFlag != timestepEnum::CourantTimeStep)
					mesh.setTimestep(
							std::min(mesh.timestep(), mesh.timestepSource()));
				mesh.setTimestep(std::max(mesh.timestep(), minTime));
				mesh.timestepSourceRef() = veryBig;

#ifdef MPI_VERSION
				/*Calculate parallel timestep*/
				{
					const MPIHandler::mpi_scalar timestep_arr[1] {
							mesh.timestep() };

					MPIHandler::mpi_scalar * nodesTimes = nullptr;
					if (parallelism.isRoot())
						nodesTimes =
								new MPIHandler::mpi_scalar[parallelism.mpi_size];

					MPI_Gather(timestep_arr, 1, MPI_DOUBLE, nodesTimes, 1,
					MPI_DOUBLE, parallelism.root, MPI_COMM_WORLD);

					MPIHandler::mpi_scalar minOfAllTimes[1] { 0 };
					if (parallelism.isRoot())
					{
						std::valarray<scalar> nodesTimes_min(
								parallelism.mpi_size);

						for (std::size_t i = 0; i < nodesTimes_min.size(); ++i)
							nodesTimes_min[i] = nodesTimes[i];

						minOfAllTimes[0] = nodesTimes_min.min();

						delete[] nodesTimes;
					}
					MPI_Bcast(minOfAllTimes, 1, MPI_DOUBLE, parallelism.root,
					MPI_COMM_WORLD);

					mesh.setTimestep(minOfAllTimes[0]);
				}
#endif
			}

			/*Write output.*/
			if ((nsteps == (nouts * frequencyOfOutput))
					|| (Time == timeOfCalculation))
			{
				std::cout << "Time = " << Time << std::endl;
				std::cout << "Time-step = " << mesh.timestep() << std::endl;
				std::cout << "Step number = " << nsteps << std::endl;

				structForOutput outputData(parallelism, mesh,
						numberOfComponents);
				if (parallelism.isRoot())
					outputData.setSizes();

				outputData.collectParallelData(*gasPhase, gasPhase->tNu,
						sonicSpeed);
				try
				{
					if (parallelism.isRoot())
						output::dataOutput(outputData, nouts, Time);
				} catch (const std::exception & exception)
				{
					std::clog << exception.what() << std::endl;

#ifndef MPI_VERSION
					inputString();
#endif

					std::exit(EXIT_FAILURE);
				}

				nouts++;
			}
			if (((nsteps == (noutsW * frequencyOfOutputWidth))
					|| (Time == timeOfCalculation))
					&& (dimensionsFlag == dimensionsEnum::task1D))
			{
				std::cout << "Time = " << Time << std::endl;
				std::cout << "Time-step = " << mesh.timestep() << std::endl;
				std::cout << "Step number = " << nsteps << std::endl;
				structForOutput outputData(parallelism, mesh,
						numberOfComponents);
				if (parallelism.isRoot())
					outputData.setSizes();

				outputData.collectParallelData(*gasPhase, gasPhase->tNu,
						sonicSpeed);
				if (parallelism.isRoot())
					output::mixedZoneWidth1D(outputData, Time);
				noutsW++;
			}

			/*Time correction for last time-step.*/
			if ((Time + mesh.timestep()) > timeOfCalculation)
				mesh.setTimestep(timeOfCalculation - Time);

#ifdef MPI_VERSION
			(MPI_Barrier(MPI_COMM_WORLD));
#endif
		} while (Time < timeOfCalculation);
		std::cout << "End" << std::endl;
#ifdef MPI_VERSION
		MPI_Finalize();
#endif

		std::ofstream outputTimeExecutionFile { "timeOfCalculation.tsv" };
		outputTimeExecutionFile.precision(ioPrecision);
		outputTimeExecutionFile << "Time of execution" << std::endl
				<< std::chrono::duration_cast<std::chrono::milliseconds>(
						std::chrono::high_resolution_clock::now() - startTime).count()
				<< std::endl << "Number of steps" << std::endl << nsteps
				<< std::endl;
		outputTimeExecutionFile << "Time of TVD stage work:" << std::endl
				<< timeForTVD << std::endl;
		outputTimeExecutionFile << "Time of Hancock stage work:" << std::endl
				<< timeForHancock << std::endl;
		outputTimeExecutionFile << "Time of Riemann solver work:" << std::endl
				<< timeForFlowCalculation << std::endl;
		outputTimeExecutionFile << "Time of time integration work:" << std::endl
				<< timeForTimeIntegration << std::endl;
		outputTimeExecutionFile << "Time of Diffusion solver work:" << std::endl
				<< timeForDiffusion << std::endl;
		outputTimeExecutionFile.close();

	} catch (const std::exception & exception)
	{
		std::clog << exception.what() << std::endl;

#ifndef MPI_VERSION
		inputString();
#endif

		std::exit(EXIT_FAILURE);
	}

	return EXIT_SUCCESS;
}
