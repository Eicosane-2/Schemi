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
#include "inputString.hpp"
#include "MPIHandler.hpp"
#include "output.hpp"
#include "phaseInitialization.hpp"
#include "scalar.hpp"
#include "secondOrderStepSolver.hpp"
#include "secondOrderStepSolverRK.hpp"
#include "structForOutput.hpp"
#include "thirdOrderStepSolver.hpp"
#include "thirdOrderStepSolverCada.hpp"
#include "ThomasSolver.hpp"
#include "typeOfSolverEnum.hpp"
#include "zone.hpp"

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
					errors::initialisationError);
#endif

		vector systemSize;
		scalar timeOfCalculation, minTime, Courant, Time { 0 };
		std::pair<scalar, scalar> timestepCoeffs { 0.25, 0.5 };
		std::size_t numberOfCells_x, numberOfCells_y, numberOfCells_z,
				frequencyOfOutput, frequencyOfOutputWidth, numberOfComponents,
				numberOfIterations, nsteps { 0 }, nouts { 0 }, noutsW { 0 };
		std::string universalGasConstant, equationOfState,
				typeOfTVDLimeterString, diffusionONString, turbulenceONString,
				matrixSolverString, gravitationONString, sourceTypeString,
				flowSolwerString, dimensionsOfTask, sourceTimeFlagString,
				linearFlagString, thirdOrderString, readFromOutput,
				mixedZoneWidthCalString;
		bool diffusionFlag, gravitationFlag, linearFlag, mixedZoneWidthCalcFlag;
		typeOfSolverEnum order;
		vector g { 0 }, gDelta { 0 };
		dimensions dimensionsFlag;
		timestep sourceTimeFlag;
		std::pair<bool, scalar> constTimeStep { false, 0. };

		std::string skipBuffer;

		/*Reading gravitational acceleration.*/
		{
			std::ifstream gFile { "./set/g.txt" };
			if (gFile.is_open())
				std::cout << "./set/g.txt is opened." << std::endl;
			else
				[[unlikely]]
				throw std::ifstream::failure("./set/g.txt not found.");

			gFile >> skipBuffer >> gravitationONString >> skipBuffer
					>> std::get<0>(g.r()) >> std::get<1>(g.r())
					>> std::get<2>(g.r()) >> skipBuffer
					>> std::get<0>(gDelta.r()) >> std::get<1>(gDelta.r())
					>> std::get<2>(gDelta.r());

			gFile.close();
		}

		/*Reading initial and boundary conditions.*/
		{
			std::ifstream mainParametersFile { "./set/main.txt" };
			if (mainParametersFile.is_open())
				std::cout << "./set/main.txt is opened." << std::endl;
			else
				[[unlikely]]
				throw std::ifstream::failure("./set/main.txt not found.");

			mainParametersFile >> skipBuffer

			>> std::get<0>(systemSize.r()) >> std::get<1>(systemSize.r())
					>> std::get<2>(systemSize.r())

					>> skipBuffer >> numberOfCells_x >> numberOfCells_y
					>> numberOfCells_z

					>> skipBuffer >> timeOfCalculation

					>> skipBuffer >> Courant

					>> skipBuffer >> frequencyOfOutput

					>> skipBuffer >> minTime

					>> skipBuffer >> numberOfComponents

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

					>> skipBuffer >> readFromOutput

					>> skipBuffer >> mixedZoneWidthCalString;

			mainParametersFile.close();
		}

		/*Setting flags.*/
		std::map<std::string, bool> reconstructionType;
		reconstructionType.insert( { "linear", true });
		reconstructionType.insert( { "star", false });
		try
		{
			linearFlag = reconstructionType.at(linearFlagString);
		} catch (const std::out_of_range&)
		{
			throw exception("Unknown reconstruction flag.",
					errors::initialisationError);
		}

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
						errors::initialisationError);

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

			std::map<std::string, timestep> sourceTimeType;
			sourceTimeType.insert( { "Courant", timestep::CourantTimeStep });
			sourceTimeType.insert( { "Courant-Source",
					timestep::CourantAndSourceTimeStep });
			sourceTimeType.insert( { "Courant-Source-Diffusion",
					timestep::CourantAndSourceAndDiffusionTimeStep });
			try
			{
				sourceTimeFlag = sourceTimeType.at(sourceTimeFlagType);
			} catch (const std::out_of_range&)
			{
				throw exception("Unknown source time-step flag.",
						errors::initialisationError);
			};
		}

		try
		{
			gravitationFlag = onOffMap.at(gravitationONString);
		} catch (const std::out_of_range&)
		{
			throw exception("Unknown gravitation flag.",
					errors::initialisationError);
		}

		try
		{
			diffusionFlag = onOffMap.at(diffusionONString);
		} catch (const std::out_of_range&)
		{
			throw exception("Unknown diffusion flag.",
					errors::initialisationError);
		}

		std::unique_ptr<abstractLimiter> limiter(
				abstractLimiter::createLimiter(typeOfTVDLimeterString));

		auto [msolver, msolverEnthFl] =
				abstractMatrixSolver::createMatrixSolver(matrixSolverString,
						dimensionsOfTask, numberOfIterations);

		std::unique_ptr<abstractFlowSolver> fsolver(
				abstractFlowSolver::createFlowSolver(flowSolwerString,
						parallelism));

		std::map<std::string, dimensions> dimensionsMap;
		dimensionsMap.insert( { "1D", dimensions::task1D });
		dimensionsMap.insert( { "2D", dimensions::task2D });
		dimensionsMap.insert( { "3D", dimensions::task3D });
		try
		{
			dimensionsFlag = dimensionsMap.at(dimensionsOfTask);
		} catch (const std::out_of_range&)
		{
			throw exception("Unknown dimensions flag.",
					errors::initialisationError);
		}

		std::map<std::string, typeOfSolverEnum> orderType;
		orderType.insert( { "ThirdOrder", typeOfSolverEnum::ThirdOrder });
		orderType.insert(
				{ "ThirdOrderCada", typeOfSolverEnum::ThirdOrderCada });
		orderType.insert( { "SecondOrder", typeOfSolverEnum::SecondOrder });
		orderType.insert( { "SecondOrderRK", typeOfSolverEnum::SecondOrderRK });
		try
		{
			order = orderType.at(thirdOrderString);
		} catch (const std::out_of_range&)
		{
			throw exception("Unknown gas dynamics approximation order.",
					errors::initialisationError);
		}

		try
		{
			mixedZoneWidthCalcFlag = onOffMap.at(mixedZoneWidthCalString);
		} catch (const std::out_of_range&)
		{
			throw exception("Unknown flag for mixed zone width calculation.",
					errors::initialisationError);
		}

		std::pair<std::size_t, std::string> readDataPoint;
		if ((readFromOutput == "no") || (readFromOutput == "initialisation"))
		{
			readDataPoint = { 0, readFromOutput };

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
				if ((dimensionsFlag == dimensions::task1D)
						&& mixedZoneWidthCalcFlag)
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
			readDataPoint = { std::stoul(readFromOutput), std::string(
					"fromTimePoint") };

			if (parallelism.isRoot())
			{
				const std::string resultFileName { "./result" };
				if (!std::filesystem::exists(resultFileName))
					throw exception(
							std::string("<<result>> directory doesn't exist."),
							errors::systemError);

				std::string timeFileName("./result/Time.tsv");
				std::ifstream timeFile(timeFileName);
				if (!timeFile.is_open())
					throw std::ifstream::failure(
							std::string("Couldn't open Time.tsv"));
				timeFile.precision(ioPrecision);

				std::string timeString;
				std::size_t lineNumber { 0 };
				std::size_t lineNumberEnd { 0 };
				bool isLastOutput { false };
				while (std::getline(timeFile, skipBuffer))
				{
					lineNumberEnd++;

					if ((lineNumberEnd - 1) == readDataPoint.first)
					{
						lineNumber = lineNumberEnd;

						timeString = skipBuffer;
					}
				}

				if (lineNumber == 0)
					throw exception(
							std::string(
									"Appropriate line wasn't found in Time.tsv"),
							errors::systemError);

				if (lineNumber == 1)
					throw exception(
							std::string("Time.tsv contains only one line."),
							errors::systemError);

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
							errors::initialisationError);

				if (nouts != readDataPoint.first)
					throw exception(
							std::string(
									"Number of output does not concur with <<readDataPoint>>."),
							errors::systemError);

				std::vector<std::string> savedStrings(lineNumber);
				timeFile.clear();
				timeFile.seekg(0, timeFile.beg);
				for (std::size_t str_i = 0; str_i < lineNumber; ++str_i)
					std::getline(timeFile, savedStrings[str_i]);

				timeFile.close();

				/*Recreation of Time.tsv.*/
				std::ofstream timeFileNew(timeFileName);
				if (!timeFileNew.is_open())
					throw std::ofstream::failure(
							std::string("Couldn't open Time.tsv"));
				timeFileNew.precision(ioPrecision);

				for (std::size_t str_i = 0; str_i < lineNumber; ++str_i)
					timeFileNew << savedStrings[str_i] << '\n';

				timeFileNew.close();

				if (isLastOutput)
				{
					try
					{
						std::ifstream numberOfStepsFile(
								"./timeOfCalculation.tsv");
						if (!numberOfStepsFile.is_open())
							throw std::ifstream::failure(
									std::string(
											"Couldn't open timeOfCalculation.tsv"));

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
						numberOfStepsFile.close();
					} catch (...)
					{
						nsteps = nouts * frequencyOfOutput;
					}
				}
				else
					nsteps = nouts * frequencyOfOutput;

				if ((dimensionsFlag == dimensions::task1D)
						&& mixedZoneWidthCalcFlag)
				{
					std::string timeWidthFileName("./result/timeWidth.tsv");
					std::ifstream timeWidthFile(timeWidthFileName);
					if (!timeWidthFile.is_open())
						throw std::ifstream::failure(
								std::string("Couldn't open timeWidth.tsv"));

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
						throw std::ofstream::failure(
								std::string("Couldn't open timeWidth.tsv"));
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
			case dimensions::task3D:
				commonConditions[2] = boundaryConditionType::calculated;
				commonConditions[5] = boundaryConditionType::calculated;
				[[fallthrough]];
			case dimensions::task2D:
				commonConditions[3] = boundaryConditionType::calculated;
				commonConditions[4] = boundaryConditionType::calculated;
				[[fallthrough]];
			case dimensions::task1D:
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
			case dimensions::task3D:
				meshObjPointer->threeDParallelepiped(parallelNodeSystemSize,
						numberOfCells_x, numberOfCells_y, numberOfCells_z,
						commonConditions);
				break;
			case dimensions::task2D:
				meshObjPointer->twoDParallelepiped(parallelNodeSystemSize,
						numberOfCells_x, numberOfCells_y, commonConditions);
				break;
			case dimensions::task1D:
			default:
				meshObjPointer->oneDParallelepiped(parallelNodeSystemSize,
						numberOfCells_x, commonConditions);
				break;
			}
		}

		auto & mesh_ = *meshObjPointer;

		parallelism.initialiseBuffersSize(mesh_);
		parallelism.initialiseParalleMeshData(mesh_);

		/*Check number of components.*/
		if (numberOfComponents > 9)
			throw exception("More than 9 components.",
					errors::initialisationError);

		const auto numberOfZones = zone::zonesArray();

		/*Creating fields.*/
		auto [gasPhase, enthalpyFlowFlag, molMassDiffusionFlag] =
				phaseInitialization(numberOfComponents, numberOfZones, mesh_,
						commonConditions, parallelism, turbulenceONString,
						sourceTypeString, universalGasConstant, equationOfState,
						readDataPoint);

		volumeField<scalar> sonicSpeed { mesh_, 0 };
		sonicSpeed.r() = std::sqrt(
				gasPhase->phaseThermodynamics->sqSonicSpeed(
						gasPhase->concentration.p, gasPhase->density[0](),
						gasPhase->internalEnergy(), gasPhase->pressure()));

		/*Calculate effective length for time-step calculation and set first time-step.*/
		volumeField<scalar> minimalLengthScale { mesh_, 0 };
		for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		{
			const auto edge1 = (mesh_.cells()[i].rX00()
					- mesh_.cells()[i].r000()).mag();
			const auto edge2 = (mesh_.cells()[i].r0Y0()
					- mesh_.cells()[i].r000()).mag();
			const auto edge3 = (mesh_.cells()[i].r00Z()
					- mesh_.cells()[i].r000()).mag();

			const auto minEdge = std::min(std::min(edge1, edge2), edge3);

			minimalLengthScale.r()[i] = 1. / minEdge;
		}

		/*Set small initial step for safe source integration.*/
		if (!constTimeStep.first)
		{
			mesh_.setTimestep(minTime);
			mesh_.setTimestep(std::min(mesh_.timestep(), timeOfCalculation));
		}
		else
			mesh_.setTimestep(constTimeStep.second);

		chemicalReactions chemReactionFlag;
		{
			std::ifstream chem { "./set/chemicalKinetics.txt" };

			if (chem.is_open())
				std::cout << "./set/chemicalKinetics.txt is opened."
						<< std::endl;
			else
				[[unlikely]]
				throw std::ifstream::failure(
						"./set/chemicalKinetics.txt not found.");

			std::string reactionName;

			chem >> skipBuffer >> reactionName;

			std::map<std::string, chemicalReactions> chemicalReactionsMap;
			chemicalReactionsMap.insert(
					{ "no", chemicalReactions::noReaction });
			chemicalReactionsMap.insert( { "Cl2",
					chemicalReactions::Cl2Dissociation });
			chemicalReactionsMap.insert( { "Cl2H2",
					chemicalReactions::Cl2H2Dissociation });
			chemicalReactionsMap.insert( { "H2Cl2Combustion",
					chemicalReactions::H2Cl2Combustion });
			chemicalReactionsMap.insert( { "NO2Disproportionation",
					chemicalReactions::NO2Disproportionation });
			chemicalReactionsMap.insert( { "H2O2Combustion",
					chemicalReactions::H2O2Combustion });
			chemicalReactionsMap.insert( { "Rober", chemicalReactions::Rober });
			try
			{
				chemReactionFlag = chemicalReactionsMap.at(reactionName);
			} catch (const std::out_of_range&)
			{
				throw exception("Unknown chemical reaction model.",
						errors::initialisationError);
			}

			chem.close();
		}

		std::unique_ptr<chemicalKinetics::abstractChemicalKinetics> chmk(
				chemicalKinetics::abstractChemicalKinetics::createChemicalKinetics(
						*gasPhase, chemReactionFlag, minTime));

		/*Write initial conditions.*/
		{
			if ((!readDataPoint.first)
					&& (readDataPoint.second != "fromTimePoint"))
			{
				structForOutput outputData(parallelism, mesh_,
						numberOfComponents);
				if (parallelism.isRoot())
					outputData.setSizes();

				outputData.collectParallelData(*gasPhase, gasPhase->tNu,
						sonicSpeed);
				if (parallelism.isRoot())
				{
					output::dataOutput(outputData, nouts, Time);

					if ((dimensionsFlag == dimensions::task1D)
							&& mixedZoneWidthCalcFlag)
						output::mixedZoneWidth1D(outputData, Time);
				}
			}
			nouts++, noutsW++;
		}
		if (readDataPoint.second == "initialisation")
			std::exit(EXIT_SUCCESS);

		boundaryConditionValue boundaryConditionValueCalc(*gasPhase->turbulence,
				*gasPhase, *(gasPhase->phaseThermodynamics), parallelism);

		skipBuffer.clear();

		parallelism.correctBoundaryValues(*gasPhase);

		std::unique_ptr<abstractStepSolver> stepSolver;

		switch (order)
		{
		case typeOfSolverEnum::SecondOrder:
			stepSolver = std::make_unique<secondOrderStepSolver>(*gasPhase,
					*limiter, *fsolver, gravitationFlag, g,
					boundaryConditionValueCalc, timeForTVD, timeForHancock,
					timeForFlowCalculation, timeForTimeIntegration, parallelism,
					diffusionFlag, *msolver, *msolverEnthFl, timestepCoeffs,
					timeForDiffusion, commonConditions, enthalpyFlowFlag,
					linearFlag, boundaryConditionValueCalc, minimalLengthScale,
					sourceTimeFlag, molMassDiffusionFlag, *chmk);
			break;
		case typeOfSolverEnum::SecondOrderRK:
			stepSolver = std::make_unique<secondOrderStepSolverRK>(*gasPhase,
					*limiter, *fsolver, gravitationFlag, g,
					boundaryConditionValueCalc, timeForTVD, timeForHancock,
					timeForFlowCalculation, timeForTimeIntegration, parallelism,
					diffusionFlag, *msolver, *msolverEnthFl, timestepCoeffs,
					timeForDiffusion, commonConditions, enthalpyFlowFlag,
					linearFlag, boundaryConditionValueCalc, minimalLengthScale,
					sourceTimeFlag, molMassDiffusionFlag, *chmk);
			break;
		case typeOfSolverEnum::ThirdOrder:
			stepSolver = std::make_unique<thirdOrderStepSolver>(*gasPhase,
					*limiter, *fsolver, gravitationFlag, g,
					boundaryConditionValueCalc, timeForTVD, timeForHancock,
					timeForFlowCalculation, timeForTimeIntegration, parallelism,
					diffusionFlag, *msolver, *msolverEnthFl, timestepCoeffs,
					timeForDiffusion, commonConditions, enthalpyFlowFlag,
					linearFlag, boundaryConditionValueCalc, minimalLengthScale,
					sourceTimeFlag, molMassDiffusionFlag, *chmk);
			break;
		case typeOfSolverEnum::ThirdOrderCada:
			stepSolver = std::make_unique<thirdOrderStepSolverCada>(*gasPhase,
					*limiter, *fsolver, gravitationFlag, g,
					boundaryConditionValueCalc, timeForTVD, timeForHancock,
					timeForFlowCalculation, timeForTimeIntegration, parallelism,
					diffusionFlag, *msolver, *msolverEnthFl, timestepCoeffs,
					timeForDiffusion, commonConditions, enthalpyFlowFlag,
					linearFlag, boundaryConditionValueCalc, minimalLengthScale,
					sourceTimeFlag, molMassDiffusionFlag, *chmk);
			break;
		[[unlikely]] default:
			throw exception("Unknown type of approximation order.",
					errors::initialisationError);
			break;
		}

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
				std::clog << "Time-step " << mesh_.timestep() << '.'
						<< std::endl;
#ifndef MPI_VERSION
				inputString();
				std::clog << "Emergency output." << std::endl;
				{
					structForOutput outputData(parallelism, mesh_,
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
				g += gDelta * mesh_.timestep();

			Time += mesh_.timestep();

			/*Calculation of new time-step through speed of sound.*/
			if (!constTimeStep.first)
			{
				sonicSpeed.r() = std::sqrt(
						gasPhase->phaseThermodynamics->sqSonicSpeed(
								gasPhase->concentration.p,
								gasPhase->density[0](),
								gasPhase->internalEnergy(),
								gasPhase->pressure()));

				std::valarray<scalar> signalSpeed(gasPhase->velocity().size());

				for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
					signalSpeed[i] = gasPhase->velocity()[i].mag()
							+ sonicSpeed()[i];

				signalSpeed *= minimalLengthScale();

				const scalar maxSignalSpeed { signalSpeed.max() };

				mesh_.setTimestep(Courant / maxSignalSpeed);
				if (sourceTimeFlag != timestep::CourantTimeStep)
					mesh_.setTimestep(
							std::min(mesh_.timestep(), mesh_.timestepSource()));
				mesh_.setTimestep(std::max(mesh_.timestep(), minTime));
				mesh_.timestepSourceRef() = veryBig;

#ifdef MPI_VERSION
				/*Calculate parallel timestep*/
				{
					const MPIHandler::mpi_scalar timestep_arr[1] {
							mesh_.timestep() };

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

					mesh_.setTimestep(minOfAllTimes[0]);
				}
#endif
			}

			/*Write output.*/
			if ((nsteps == (nouts * frequencyOfOutput))
					|| (Time == timeOfCalculation))
			{
				std::cout << "Time = " << Time << std::endl;
				std::cout << "Time-step = " << mesh_.timestep() << std::endl;
				std::cout << "Step number = " << nsteps << std::endl;

				structForOutput outputData(parallelism, mesh_,
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
					&& (dimensionsFlag == dimensions::task1D)
					&& mixedZoneWidthCalcFlag)
			{
				std::cout << "Time = " << Time << std::endl;
				std::cout << "Time-step = " << mesh_.timestep() << std::endl;
				std::cout << "Step number = " << nsteps << std::endl;
				structForOutput outputData(parallelism, mesh_,
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
			if ((Time + mesh_.timestep()) > timeOfCalculation)
				mesh_.setTimestep(timeOfCalculation - Time);

#ifdef MPI_VERSION
			(MPI_Barrier(MPI_COMM_WORLD));
#endif
		} while (Time < timeOfCalculation);
		std::cout << "End" << std::endl;
#ifdef MPI_VERSION
		MPI_Finalize();
#endif

		std::ofstream outputTimeExecutionFile { "./timeOfCalculation.tsv" };
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
