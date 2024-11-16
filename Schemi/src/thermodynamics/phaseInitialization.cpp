/*
 * phaseInitialization.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "phaseInitialization.hpp"

#include <iostream>
#include <filesystem>

#include "gasModelEnum.hpp"
#include "boundaryConditionFromString.hpp"
#include "mixtureIdeal.hpp"
#include "mixtureKataokaVanDerWaals.hpp"
#include "mixtureRedlichKwong.hpp"
#include "mixtureStiffened.hpp"
#include "mixtureVanDerWaals.hpp"
#include "fieldProducts.hpp"

std::tuple<std::unique_ptr<schemi::homogeneousPhase<schemi::cubicCell>>,
		schemi::enthalpyFlow, bool> schemi::phaseInitialization(
		std::size_t numberOfComponents,
		const std::vector<std::unique_ptr<zone>> & numberOfZones,
		const mesh & meshReference,
		const std::vector<boundaryConditionType> & commonConditions,
		const MPIHandler & parallelism, const std::string & turbulenceONString,
		const std::string & sourceTypeString,
		const std::string & universalGasConstant,
		const std::string & equationOfState,
		const std::pair<std::size_t, std::string> & readDataPoint)
{
	std::unique_ptr<homogeneousPhase<cubicCell>> phase { nullptr };
	enthalpyFlow enthalpyFlowFlag;
	bool molMassDiffusionFlag;

	std::string skipBuffer;

	std::unique_ptr<abstractTurbulenceGen> turbulenceSources(
			abstractTurbulenceGen::createTurbulenceModel(meshReference,
					turbulenceONString, sourceTypeString));

	transportCoefficients<cubicCell> tCoeffsPhase { meshReference };

	bunchOfFields<cubicCell> cellFields { meshReference, numberOfComponents };
	std::unique_ptr<abstractMixtureThermodynamics> mixture;

	/*Set boundary condition for transport coefficients.*/
	auto bndConCom = commonConditions;

	std::replace(bndConCom.begin(), bndConCom.end(),
			boundaryConditionType::calculated,
			boundaryConditionType::calculatedTurbulentViscosity);
	tCoeffsPhase.tNu = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });

	std::ifstream transportCoefficientsFile { "./set/transportCoefficients.txt" };
	if (transportCoefficientsFile.is_open())
		std::cout << "./set/transportCoefficients.txt is opened." << std::endl;
	else
		throw std::ifstream::failure(
				"./set/transportCoefficients.txt not found.");

	std::string transportModelTypeString;
	scalar cNu, cD, cKappa;
	std::string implicitEnthalpyFlowFlagString, MolMassDiffusionFlagString;
	transportCoefficientsFile >> skipBuffer >> transportModelTypeString
			>> skipBuffer >> cNu >> skipBuffer >> cD >> skipBuffer >> cKappa
			>> skipBuffer >> implicitEnthalpyFlowFlagString >> skipBuffer
			>> MolMassDiffusionFlagString;
	transportCoefficientsFile.close();

	std::map<std::string, transportModel> transportModels;
	transportModels.insert( { "constant", transportModel::constant });
	transportModels.insert( { "hardSpheres", transportModel::hardSpheres });
	transportModel transpModel;
	try
	{
		transpModel = transportModels.at(transportModelTypeString);
	} catch (const std::out_of_range&)
	{
		throw exception("Unknown flag for transport model.",
				errors::initialisationError);
	}

	std::map<std::string, enthalpyFlow> enthalpyFlowTypes;
	enthalpyFlowTypes.insert( { "implicit", enthalpyFlow::implicitSolve });
	enthalpyFlowTypes.insert( { "explicit", enthalpyFlow::explicitSolve });
	enthalpyFlowTypes.insert( { "no", enthalpyFlow::noSolve });
	try
	{
		enthalpyFlowFlag = enthalpyFlowTypes.at(implicitEnthalpyFlowFlagString);
	} catch (const std::out_of_range&)
	{
		throw exception("Unknown flag for enthalpy flow calculation.",
				errors::initialisationError);
	}

	try
	{
		molMassDiffusionFlag = onOffMap.at(MolMassDiffusionFlagString);
	} catch (const std::out_of_range&)
	{
		throw exception("Unknown flag for molar mass diffusion correction.",
				errors::initialisationError);
	}

	/*Create mixture class and boundary conditions for basic quantities: concentrations, velocity, pressure, k, epsilon, a and b.*/
	std::vector<std::array<std::vector<subPatchData<scalar>>, 6>> boundaryConditionsMatrix(
			numberOfComponents);
	for (auto & arr : boundaryConditionsMatrix)
		arr.fill(std::vector<subPatchData<scalar>>(0));

	std::vector<std::vector<std::string>> matrixOfSubstancesConditions(
			numberOfComponents,
			std::vector<std::string>(5 + numberOfZones.size()));

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		std::string substanceName { "./set/sub_" }, bufComponent(
				std::to_string(k + 1));

		substanceName.append(bufComponent);
		substanceName.append(".txt");
		std::ifstream substanceConditionsFile { substanceName };

		if (substanceConditionsFile.is_open())
			std::cout << substanceName << " is opened." << std::endl;
		else
			throw std::ifstream::failure(substanceName + " not found.");

		substanceConditionsFile >> skipBuffer
				>> matrixOfSubstancesConditions[k][0] /*M*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][1] /*Cv*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][2] /*Tcrit*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][3] /*Pcrit*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][4]; /*Molecular diameter*/

		/*Read boundary*/
		for (std::size_t boundaryI = 0;
				boundaryI < boundaryConditionsMatrix[k].size(); ++boundaryI)
		{
			substanceConditionsFile >> skipBuffer;
			std::string patchBoundary;
			substanceConditionsFile >> patchBoundary;
			if (patchBoundary != "subPatches")
			{
				const auto boundTp = boundaryConditionFromString(patchBoundary);

				if (boundTp == boundaryConditionType::fixedValueCell)
				{
					scalar fv;
					substanceConditionsFile >> fv;

					boundaryConditionsMatrix[k][boundaryI].emplace_back(boundTp,
							fv);
				}
				else
					boundaryConditionsMatrix[k][boundaryI].emplace_back(
							boundTp);
			}
			else
			{
				std::size_t subPatchNumber;

				substanceConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errors::initialisationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					substanceConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValueCell)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						substanceConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						substanceConditionsFile >> fv;

						boundaryConditionsMatrix[k][boundaryI].emplace_back(
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 }, fv);
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						substanceConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsMatrix[k][boundaryI].emplace_back(
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 });
					}
				}
			}
		}

		substanceConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errors::initialisationError);
		for (std::size_t j = 0; j < numberOfZones.size(); ++j)
			substanceConditionsFile >> matrixOfSubstancesConditions[k][5 + j]; /*Values in zones*/
		substanceConditionsFile.close();
	}

	/*Convert strings to values.*/
	std::array<std::valarray<scalar>, 4> thermodynamicalProperties;
	thermodynamicalProperties.fill(std::valarray<scalar>(numberOfComponents));

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		std::get<0>(thermodynamicalProperties)[k] = std::stod(
				matrixOfSubstancesConditions[k][0]); /*M*/
		std::get<1>(thermodynamicalProperties)[k] = std::stod(
				matrixOfSubstancesConditions[k][1]); /*Cv*/
		std::get<2>(thermodynamicalProperties)[k] = std::stod(
				matrixOfSubstancesConditions[k][2]); /*Tcrit*/
		std::get<3>(thermodynamicalProperties)[k] = std::stod(
				matrixOfSubstancesConditions[k][3]); /*Pcrit*/
	}

	for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
	{
		parallelism.correctBoundaryConditions(boundaryConditionsMatrix[k - 1]);

		cellFields.concentration.v[k] = volumeField<scalar>(meshReference,
				scalar { 0 }, boundaryConditionsMatrix[k - 1][0],
				boundaryConditionsMatrix[k - 1][1],
				boundaryConditionsMatrix[k - 1][2],
				boundaryConditionsMatrix[k - 1][3],
				boundaryConditionsMatrix[k - 1][4],
				boundaryConditionsMatrix[k - 1][5]);

		auto bndVConRho_k = boundaryConditionsMatrix[k - 1];

		for (auto & b_s : bndVConRho_k)
			for (auto & b_i : b_s)
				b_i.fixVal *= thermodynamicalProperties[0][k - 1];

		cellFields.density[k] = volumeField<scalar>(meshReference, scalar { 0 },
				bndVConRho_k[0], bndVConRho_k[1], bndVConRho_k[2],
				bndVConRho_k[3], bndVConRho_k[4], bndVConRho_k[5]);
	}
	{
		cellFields.concentration.v[0] = volumeField<scalar>(meshReference, 0);

		auto bndConRho = commonConditions;

		std::replace(bndConRho.begin(), bndConRho.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedDensity);

		cellFields.density[0] = volumeField<scalar>(meshReference, scalar { 0 },
				subPatchData<scalar> { bndConRho[0] }, subPatchData<scalar> {
						bndConRho[1] }, subPatchData<scalar> { bndConRho[2] },
				subPatchData<scalar> { bndConRho[3] }, subPatchData<scalar> {
						bndConRho[4] }, subPatchData<scalar> { bndConRho[5] });
	}

	scalar R { 0 }, hPlanck { 0 };
	if (universalGasConstant == "SI")
	{
		R = R_SI;
		hPlanck = PlanckConstant_SI;
	}
	else if (universalGasConstant == "1")
	{
		R = 1;
		hPlanck = PlanckConstant_SI / R_SI;
	}
	else if (universalGasConstant == "CGS")
	{
		R = R_SI * 1E7;
		hPlanck = PlanckConstant_SI * 1E7;
	}
	else
		throw exception("Unknown type of universal gas constant.",
				errors::initialisationError);

	{
		std::map<std::string, gasModel> gasModels;
		gasModels.insert( { "vanDerWaals", gasModel::vanDerWaals });
		gasModels.insert( { "ideal", gasModel::ideal });
		gasModels.insert( { "RedlichKwong", gasModel::RedlichKwong });
		gasModels.insert( { "stiffened", gasModel::stiffened });
		gasModels.insert( { "Kataoka", gasModel::KataokaVanDerWaals });
		gasModel gasModelFlag;
		try
		{
			gasModelFlag = gasModels.at(equationOfState);
		} catch (const std::out_of_range&)
		{
			throw exception("Unknown type of equation of state.",
					errors::initialisationError);
		}

		switch (gasModelFlag)
		{
		case gasModel::vanDerWaals:
			mixture = std::make_unique<mixtureVanDerWaals>(R, hPlanck,
					std::get<0>(thermodynamicalProperties),
					std::get<1>(thermodynamicalProperties),
					std::get<2>(thermodynamicalProperties),
					std::get<3>(thermodynamicalProperties));
			break;
		case gasModel::RedlichKwong:
			mixture = std::make_unique<mixtureRedlichKwong>(R, hPlanck,
					std::get<0>(thermodynamicalProperties),
					std::get<1>(thermodynamicalProperties),
					std::get<2>(thermodynamicalProperties),
					std::get<3>(thermodynamicalProperties));
			break;
		case gasModel::stiffened:
		{
			const std::string stiffenedCoeffsName {
					"./set/stiffenedFluidData.txt" };

			std::ifstream fluidConditionsFile { stiffenedCoeffsName };

			if (fluidConditionsFile.is_open())
				std::cout << stiffenedCoeffsName << " is opened." << std::endl;
			else
				throw std::ifstream::failure(
						stiffenedCoeffsName + " not found.");

			std::valarray<scalar> p0Data(numberOfComponents), gammaData(
					numberOfComponents);

			for (std::size_t k = 0; k < numberOfComponents; ++k)
			{
				if (fluidConditionsFile.eof())
					throw std::ifstream::failure(
							stiffenedCoeffsName + ". Unexpected end of file.");

				fluidConditionsFile >> skipBuffer >> p0Data[k] >> gammaData[k];
			}

			mixture = std::make_unique<mixtureStiffened>(R, hPlanck,
					std::get<0>(thermodynamicalProperties),
					std::get<1>(thermodynamicalProperties), p0Data, gammaData);
		}
			break;
		case gasModel::KataokaVanDerWaals:
		{
			const std::string KataokaCoeffsName {
					"./set/KataokaVanDerWaalsFluidData.txt" };

			std::ifstream fluidConditionsFile { KataokaCoeffsName };

			if (fluidConditionsFile.is_open())
				std::cout << KataokaCoeffsName << " is opened." << std::endl;
			else
				throw std::ifstream::failure(KataokaCoeffsName + " not found.");

			std::valarray<scalar> epsilonLJData(numberOfComponents),
					sigmaLJData(numberOfComponents);

			for (std::size_t k = 0; k < numberOfComponents; ++k)
			{
				if (fluidConditionsFile.eof())
					throw std::ifstream::failure(
							KataokaCoeffsName + ". Unexpected end of file.");

				fluidConditionsFile >> skipBuffer >> epsilonLJData[k]
						>> sigmaLJData[k];
			}

			std::string bCalcTypeStr;
			scalar bCoeff;

			fluidConditionsFile >> skipBuffer >> bCalcTypeStr >> bCoeff;

			std::pair<bool, scalar> bCalc;

			if (bCalcTypeStr == "critical")
				bCalc = { true, bCoeff };
			else if (bCalcTypeStr == "Kataoka")
				bCalc = { false, bCoeff };
			else
				throw exception(
						"Wrong type of Kataoka-van der Waals b coefficient calculation",
						errors::initialisationError);

			mixture = std::make_unique<mixtureKataokaVanDerWaals>(R, hPlanck,
					std::get<0>(thermodynamicalProperties),
					std::get<1>(thermodynamicalProperties),
					std::get<2>(thermodynamicalProperties),
					std::get<3>(thermodynamicalProperties), epsilonLJData,
					sigmaLJData, bCalc);
		}
			break;
		case gasModel::ideal:
		default:
			mixture = std::make_unique<mixtureIdeal>(R, hPlanck,
					std::get<0>(thermodynamicalProperties),
					std::get<1>(thermodynamicalProperties));
			break;
		}
	}

	std::unique_ptr<abstractTransportModel> transportModel(
			abstractTransportModel::createTransportModel(
					matrixOfSubstancesConditions, cNu, cD, cKappa,
					transpModel));

	phase = std::make_unique<homogeneousPhase<cubicCell>>(cellFields,
			tCoeffsPhase, mixture, turbulenceSources, transportModel);

	std::array<std::vector<subPatchData<vector>>, 6> boundaryConditionsVelocity;
	boundaryConditionsVelocity.fill(std::vector<subPatchData<vector>>(0));
	std::vector<std::string> velocityConditions(3 * numberOfZones.size());
	{
		std::ifstream velocityConditionsFile { "./set/velocity.txt" };

		if (velocityConditionsFile.is_open())
			std::cout << "./set/velocity.txt is opened." << std::endl;
		else
			throw std::ifstream::failure("./set/velocity.txt not found.");

		/*Read boundary*/
		for (std::size_t boundaryI = 0;
				boundaryI < boundaryConditionsVelocity.size(); ++boundaryI)
		{
			velocityConditionsFile >> skipBuffer;
			std::string patchBoundary;
			velocityConditionsFile >> patchBoundary;
			if (patchBoundary != "subPatches")
			{
				const auto boundTp = boundaryConditionFromString(patchBoundary);

				if (boundTp == boundaryConditionType::fixedValueCell)
				{
					scalar fv1, fv2, fv3;
					velocityConditionsFile >> fv1 >> fv2 >> fv3;

					boundaryConditionsVelocity[boundaryI].emplace_back(boundTp,
							vector { fv1, fv2, fv3 });
				}
				else
					boundaryConditionsVelocity[boundaryI].emplace_back(boundTp);
			}
			else
			{
				std::size_t subPatchNumber;

				velocityConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errors::initialisationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					velocityConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValueCell)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						velocityConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv1, fv2, fv3;
						velocityConditionsFile >> fv1 >> fv2 >> fv3;

						boundaryConditionsVelocity[boundaryI].emplace_back(
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 }, vector { fv1, fv2, fv3 });
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						velocityConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsVelocity[boundaryI].emplace_back(
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 });
					}
				}
			}
		}

		velocityConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errors::initialisationError);
		for (std::size_t i = 0; i < numberOfZones.size(); ++i)
			velocityConditionsFile >> velocityConditions[3 * i]
					>> velocityConditions[1 + 3 * i]
					>> velocityConditions[2 + 3 * i];

		velocityConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsVelocity);

		phase->velocity = volumeField<vector>(meshReference, vector(0),
				std::get<0>(boundaryConditionsVelocity),
				std::get<1>(boundaryConditionsVelocity),
				std::get<2>(boundaryConditionsVelocity),
				std::get<3>(boundaryConditionsVelocity),
				std::get<4>(boundaryConditionsVelocity),
				std::get<5>(boundaryConditionsVelocity));
	}

	std::array<std::vector<subPatchData<scalar>>, 6> boundaryConditionsPressure;
	boundaryConditionsPressure.fill(std::vector<subPatchData<scalar>>(0));
	std::vector<std::string> pressureConditions(numberOfZones.size());
	{
		std::ifstream pressureConditionsFile { "./set/pressure.txt" };

		if (pressureConditionsFile.is_open())
			std::cout << "./set/pressure.txt is opened." << std::endl;
		else
			throw std::ifstream::failure("./set/pressure.txt not found.");

		/*Read boundary*/
		for (std::size_t boundaryI = 0;
				boundaryI < boundaryConditionsPressure.size(); ++boundaryI)
		{
			pressureConditionsFile >> skipBuffer;
			std::string patchBoundary;
			pressureConditionsFile >> patchBoundary;
			if (patchBoundary != "subPatches")
			{
				const auto boundTp = boundaryConditionFromString(patchBoundary);

				if (boundTp == boundaryConditionType::fixedValueCell)
				{
					scalar fv;
					pressureConditionsFile >> fv;

					boundaryConditionsPressure[boundaryI].emplace_back(boundTp,
							fv);
				}
				else
					boundaryConditionsPressure[boundaryI].emplace_back(boundTp);
			}
			else
			{
				std::size_t subPatchNumber;

				pressureConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errors::initialisationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					pressureConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValueCell)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						pressureConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						pressureConditionsFile >> fv;

						boundaryConditionsPressure[boundaryI].emplace_back(
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 }, fv);
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						pressureConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsPressure[boundaryI].emplace_back(
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 });
					}
				}
			}
		}

		pressureConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errors::initialisationError);
		for (std::size_t i = 0; i < numberOfZones.size(); ++i)
			pressureConditionsFile >> pressureConditions[i];

		pressureConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsPressure);

		phase->pressure = volumeField<scalar>(meshReference, scalar { 0 },
				std::get<0>(boundaryConditionsPressure),
				std::get<1>(boundaryConditionsPressure),
				std::get<2>(boundaryConditionsPressure),
				std::get<3>(boundaryConditionsPressure),
				std::get<4>(boundaryConditionsPressure),
				std::get<5>(boundaryConditionsPressure));
	}

	std::array<std::vector<subPatchData<scalar>>, 6> boundaryConditionskTurb;
	boundaryConditionskTurb.fill(std::vector<subPatchData<scalar>>(0));
	std::vector<std::string> kTurbConditions(numberOfZones.size());
	{
		std::ifstream kTurbConditionsFile { "./set/k.txt" };

		if (kTurbConditionsFile.is_open())
			std::cout << "./set/k.txt is opened." << std::endl;
		else
			throw std::ifstream::failure("./set/k.txt is opened.");

		/*Read boundary*/
		for (std::size_t boundaryI = 0;
				boundaryI < boundaryConditionskTurb.size(); ++boundaryI)
		{
			kTurbConditionsFile >> skipBuffer;
			std::string patchBoundary;
			kTurbConditionsFile >> patchBoundary;
			if (patchBoundary != "subPatches")
			{
				const auto boundTp = boundaryConditionFromString(patchBoundary);

				if (boundTp == boundaryConditionType::fixedValueCell)
				{
					scalar fv;
					kTurbConditionsFile >> fv;

					boundaryConditionskTurb[boundaryI].emplace_back(boundTp,
							fv);
				}
				else
					boundaryConditionskTurb[boundaryI].emplace_back(boundTp);
			}
			else
			{
				std::size_t subPatchNumber;

				kTurbConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errors::initialisationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					kTurbConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValueCell)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						kTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						kTurbConditionsFile >> fv;

						boundaryConditionskTurb[boundaryI].emplace_back(boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 }, fv);
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						kTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionskTurb[boundaryI].emplace_back(boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 });
					}
				}
			}
		}

		kTurbConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errors::initialisationError);
		for (std::size_t i = 0; i < numberOfZones.size(); ++i)
			kTurbConditionsFile >> kTurbConditions[i];

		kTurbConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionskTurb);

		phase->kTurb = volumeField<scalar>(meshReference, scalar { 0 },
				std::get<0>(boundaryConditionskTurb),
				std::get<1>(boundaryConditionskTurb),
				std::get<2>(boundaryConditionskTurb),
				std::get<3>(boundaryConditionskTurb),
				std::get<4>(boundaryConditionskTurb),
				std::get<5>(boundaryConditionskTurb));
	}

	std::array<std::vector<subPatchData<scalar>>, 6> boundaryConditionsepsTurb;
	boundaryConditionsepsTurb.fill(std::vector<subPatchData<scalar>>(0));
	std::vector<std::string> epsTurbConditions(numberOfZones.size());
	{
		std::ifstream epsTurbConditionsFile { "./set/epsilon.txt" };

		if (epsTurbConditionsFile.is_open())
			std::cout << "./set/epsilon.txt is opened." << std::endl;
		else
			throw std::ifstream::failure("./set/epsilon.txt not found.");

		/*Read boundary*/
		for (std::size_t boundaryI = 0;
				boundaryI < boundaryConditionsepsTurb.size(); ++boundaryI)
		{
			epsTurbConditionsFile >> skipBuffer;
			std::string patchBoundary;
			epsTurbConditionsFile >> patchBoundary;
			if (patchBoundary != "subPatches")
			{
				const auto boundTp = boundaryConditionFromString(patchBoundary);

				if (boundTp == boundaryConditionType::fixedValueCell)
				{
					scalar fv;
					epsTurbConditionsFile >> fv;

					boundaryConditionsepsTurb[boundaryI].emplace_back(boundTp,
							fv);
				}
				else
					boundaryConditionsepsTurb[boundaryI].emplace_back(boundTp);
			}
			else
			{
				std::size_t subPatchNumber;

				epsTurbConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errors::initialisationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					epsTurbConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValueCell)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						epsTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						epsTurbConditionsFile >> fv;

						boundaryConditionsepsTurb[boundaryI].emplace_back(
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 }, fv);
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						epsTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsepsTurb[boundaryI].emplace_back(
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 });
					}
				}
			}
		}

		epsTurbConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errors::initialisationError);
		for (std::size_t i = 0; i < numberOfZones.size(); ++i)
			epsTurbConditionsFile >> epsTurbConditions[i];

		epsTurbConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsepsTurb);

		phase->epsTurb = volumeField<scalar>(meshReference, scalar { 0 },
				std::get<0>(boundaryConditionsepsTurb),
				std::get<1>(boundaryConditionsepsTurb),
				std::get<2>(boundaryConditionsepsTurb),
				std::get<3>(boundaryConditionsepsTurb),
				std::get<4>(boundaryConditionsepsTurb),
				std::get<5>(boundaryConditionsepsTurb));
	}

	std::array<std::vector<subPatchData<vector>>, 6> boundaryConditionsaTurb;
	boundaryConditionsaTurb.fill(std::vector<subPatchData<vector>>(0));
	std::vector<std::string> aTurbConditions(3 * numberOfZones.size());
	{
		std::ifstream aTurbConditionsFile { "./set/a.txt" };

		if (aTurbConditionsFile.is_open())
			std::cout << "./set/a.txt is opened." << std::endl;
		else
			throw std::ifstream::failure("./set/a.txt not found.");

		/*Read boundary*/
		for (std::size_t boundaryI = 0;
				boundaryI < boundaryConditionsaTurb.size(); ++boundaryI)
		{
			aTurbConditionsFile >> skipBuffer;
			std::string patchBoundary;
			aTurbConditionsFile >> patchBoundary;
			if (patchBoundary != "subPatches")
			{
				const auto boundTp = boundaryConditionFromString(patchBoundary);

				if (boundTp == boundaryConditionType::fixedValueCell)
				{
					scalar fv1, fv2, fv3;
					aTurbConditionsFile >> fv1 >> fv2 >> fv3;

					boundaryConditionsaTurb[boundaryI].emplace_back(boundTp,
							vector { fv1, fv2, fv3 });
				}
				else
					boundaryConditionsaTurb[boundaryI].emplace_back(boundTp);
			}
			else
			{
				std::size_t subPatchNumber;

				aTurbConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errors::initialisationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					aTurbConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValueCell)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						aTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv1, fv2, fv3;
						aTurbConditionsFile >> fv1 >> fv2 >> fv3;

						boundaryConditionsaTurb[boundaryI].emplace_back(boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 },
								vector { fv1, fv2, fv3 });
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						aTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsaTurb[boundaryI].emplace_back(boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 });
					}
				}
			}
		}

		aTurbConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errors::initialisationError);
		for (std::size_t i = 0; i < numberOfZones.size(); ++i)
			aTurbConditionsFile >> aTurbConditions[3 * i]
					>> aTurbConditions[1 + 3 * i] >> aTurbConditions[2 + 3 * i];

		aTurbConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsaTurb);

		phase->aTurb = volumeField<vector>(meshReference, vector(0, 0, 0),
				std::get<0>(boundaryConditionsaTurb),
				std::get<1>(boundaryConditionsaTurb),
				std::get<2>(boundaryConditionsaTurb),
				std::get<3>(boundaryConditionsaTurb),
				std::get<4>(boundaryConditionsaTurb),
				std::get<5>(boundaryConditionsaTurb));
	}

	std::array<std::vector<subPatchData<scalar>>, 6> boundaryConditionsbTurb;
	boundaryConditionsbTurb.fill(std::vector<subPatchData<scalar>>(0));
	std::vector<std::string> bTurbConditions(numberOfZones.size());
	{
		std::ifstream bTurbConditionsFile { "./set/b.txt" };

		if (bTurbConditionsFile.is_open())
			std::cout << "./set/b.txt is opened." << std::endl;
		else
			throw std::ifstream::failure("./set/b.txt not found.");

		/*Read boundary*/
		for (std::size_t boundaryI = 0;
				boundaryI < boundaryConditionsbTurb.size(); ++boundaryI)
		{
			bTurbConditionsFile >> skipBuffer;
			std::string patchBoundary;
			bTurbConditionsFile >> patchBoundary;
			if (patchBoundary != "subPatches")
			{
				const auto boundTp = boundaryConditionFromString(patchBoundary);

				if (boundTp == boundaryConditionType::fixedValueCell)
				{
					scalar fv;
					bTurbConditionsFile >> fv;

					boundaryConditionsbTurb[boundaryI].emplace_back(boundTp,
							fv);
				}
				else
					boundaryConditionsbTurb[boundaryI].emplace_back(boundTp);
			}
			else
			{
				std::size_t subPatchNumber;

				bTurbConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errors::initialisationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					bTurbConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValueCell)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						bTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						bTurbConditionsFile >> fv;

						boundaryConditionsbTurb[boundaryI].emplace_back(boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 }, fv);
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						bTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsbTurb[boundaryI].emplace_back(boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 });
					}
				}
			}
		}

		bTurbConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errors::initialisationError);
		for (std::size_t i = 0; i < numberOfZones.size(); ++i)
			bTurbConditionsFile >> bTurbConditions[i];

		bTurbConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsbTurb);

		phase->bTurb = volumeField<scalar>(meshReference,
				phase->turbulenceSources->turbPar->minb_value,
				std::get<0>(boundaryConditionsbTurb),
				std::get<1>(boundaryConditionsbTurb),
				std::get<2>(boundaryConditionsbTurb),
				std::get<3>(boundaryConditionsbTurb),
				std::get<4>(boundaryConditionsbTurb),
				std::get<5>(boundaryConditionsbTurb));
	}

	/*Creating other fields.*/
	phase->momentum = volumeField<vector>(meshReference, vector(0),
			subPatchData<vector> { commonConditions[0] }, subPatchData<vector> {
					commonConditions[1] }, subPatchData<vector> {
					commonConditions[2] }, subPatchData<vector> {
					commonConditions[3] }, subPatchData<vector> {
					commonConditions[4] }, subPatchData<vector> {
					commonConditions[5] });

	phase->internalEnergy = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { commonConditions[0] }, subPatchData<scalar> {
					commonConditions[1] }, subPatchData<scalar> {
					commonConditions[2] }, subPatchData<scalar> {
					commonConditions[3] }, subPatchData<scalar> {
					commonConditions[4] }, subPatchData<scalar> {
					commonConditions[5] });

	phase->totalEnergy = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { commonConditions[0] }, subPatchData<scalar> {
					commonConditions[1] }, subPatchData<scalar> {
					commonConditions[2] }, subPatchData<scalar> {
					commonConditions[3] }, subPatchData<scalar> {
					commonConditions[4] }, subPatchData<scalar> {
					commonConditions[5] });

	{
		auto bndCon = commonConditions;
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedTemperature);
		phase->temperature = volumeField<scalar>(meshReference, 0,
				subPatchData<scalar> { bndCon[0] }, subPatchData<scalar> {
						bndCon[1] }, subPatchData<scalar> { bndCon[2] },
				subPatchData<scalar> { bndCon[3] }, subPatchData<scalar> {
						bndCon[4] }, subPatchData<scalar> { bndCon[5] });
	}

	phase->rhokTurb = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { commonConditions[0] }, subPatchData<scalar> {
					commonConditions[1] }, subPatchData<scalar> {
					commonConditions[2] }, subPatchData<scalar> {
					commonConditions[3] }, subPatchData<scalar> {
					commonConditions[4] }, subPatchData<scalar> {
					commonConditions[5] });

	phase->rhoepsTurb = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { commonConditions[0] }, subPatchData<scalar> {
					commonConditions[1] }, subPatchData<scalar> {
					commonConditions[2] }, subPatchData<scalar> {
					commonConditions[3] }, subPatchData<scalar> {
					commonConditions[4] }, subPatchData<scalar> {
					commonConditions[5] });

	phase->rhoaTurb = volumeField<vector>(meshReference, vector(0),
			subPatchData<vector> { commonConditions[0] }, subPatchData<vector> {
					commonConditions[1] }, subPatchData<vector> {
					commonConditions[2] }, subPatchData<vector> {
					commonConditions[3] }, subPatchData<vector> {
					commonConditions[4] }, subPatchData<vector> {
					commonConditions[5] });

	phase->rhobTurb = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { commonConditions[0] }, subPatchData<scalar> {
					commonConditions[1] }, subPatchData<scalar> {
					commonConditions[2] }, subPatchData<scalar> {
					commonConditions[3] }, subPatchData<scalar> {
					commonConditions[4] }, subPatchData<scalar> {
					commonConditions[5] });

	phase->HelmholtzEnergy = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { commonConditions[0] }, subPatchData<scalar> {
					commonConditions[1] }, subPatchData<scalar> {
					commonConditions[2] }, subPatchData<scalar> {
					commonConditions[3] }, subPatchData<scalar> {
					commonConditions[4] }, subPatchData<scalar> {
					commonConditions[5] });

	phase->entropy = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { commonConditions[0] }, subPatchData<scalar> {
					commonConditions[1] }, subPatchData<scalar> {
					commonConditions[2] }, subPatchData<scalar> {
					commonConditions[3] }, subPatchData<scalar> {
					commonConditions[4] }, subPatchData<scalar> {
					commonConditions[5] });

	if ((readDataPoint.second == "no")
			|| (readDataPoint.second == "initialisation"))
	{
		/*Set initial conditions for concentrations.*/
		for (std::size_t k = 1; k < phase->concentration.v.size(); ++k)
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				const vector & radiusOfCell { meshReference.cells()[i].rC() };
				bool zoneFounded { false };

				for (std::size_t j = 0; j < numberOfZones.size(); ++j)
				{
					zoneFounded = numberOfZones[j]->isPointInZone(radiusOfCell);

					if (zoneFounded)
					{
						phase->concentration.v[k].r()[i] = std::stod(
								matrixOfSubstancesConditions[k - 1][5 + j]);
						break;
					}
				}

				if (!zoneFounded)
					throw exception(
							"Cell " + std::to_string(i) + " of "
									+ std::to_string(k)
									+ " component concentration dropped out of all zones.",
							errors::initialisationError);
			}

		/*Set initial conditions for velocity.*/
		for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
		{
			const vector & radiusOfCell { meshReference.cells()[i].rC() };
			bool zoneFounded { false };

			for (std::size_t j = 0; j < numberOfZones.size(); ++j)
			{
				zoneFounded = numberOfZones[j]->isPointInZone(radiusOfCell);

				if (zoneFounded)
				{
					phase->velocity.r()[i] = vector(
							std::stod(velocityConditions[3 * j]),
							std::stod(velocityConditions[1 + 3 * j]),
							std::stod(velocityConditions[2 + 3 * j]));
					break;
				}
			}

			if (!zoneFounded)
				throw exception(
						"Cell " + std::to_string(i)
								+ " of velocity dropped out of all zones.",
						errors::initialisationError);
		}

		/*Set initial conditions for pressure.*/
		for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
		{
			const vector & radiusOfCell { meshReference.cells()[i].rC() };
			bool zoneFounded { false };

			for (std::size_t j = 0; j < numberOfZones.size(); ++j)
			{
				zoneFounded = numberOfZones[j]->isPointInZone(radiusOfCell);

				if (zoneFounded)
				{
					phase->pressure.r()[i] = std::stod(pressureConditions[j]);
					break;
				}
			}
			if (!zoneFounded)
				throw exception(
						"Cell " + std::to_string(i)
								+ " of pressure dropped out of all zones.",
						errors::initialisationError);
		}

		/*Set initial conditions for k.*/
		for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
		{
			const vector & radiusOfCell { meshReference.cells()[i].rC() };
			bool zoneFounded { false };

			for (std::size_t j = 0; j < numberOfZones.size(); ++j)
			{
				zoneFounded = numberOfZones[j]->isPointInZone(radiusOfCell);

				if (zoneFounded)
				{
					phase->kTurb.r()[i] = std::stod(kTurbConditions[j]);
					break;
				}
			}
			if (!zoneFounded)
				throw exception(
						"Cell " + std::to_string(i)
								+ " of k dropped out of all zones.",
						errors::initialisationError);
		}

		/*Set initial conditions for epsilon.*/
		for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
		{
			const vector & radiusOfCell { meshReference.cells()[i].rC() };
			bool zoneFounded { false };

			for (std::size_t j = 0; j < numberOfZones.size(); ++j)
			{
				zoneFounded = numberOfZones[j]->isPointInZone(radiusOfCell);

				if (zoneFounded)
				{
					phase->epsTurb.r()[i] = std::stod(epsTurbConditions[j]);
					break;
				}
			}
			if (!zoneFounded)
				throw exception(
						"Cell " + std::to_string(i)
								+ " of epsilon dropped out of all zones.",
						errors::initialisationError);
		}

		/*Set initial conditions for a.*/
		for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
		{
			const vector & radiusOfCell { meshReference.cells()[i].rC() };
			bool zoneFounded { false };

			for (std::size_t j = 0; j < numberOfZones.size(); ++j)
			{
				zoneFounded = numberOfZones[j]->isPointInZone(radiusOfCell);

				if (zoneFounded)
				{
					phase->aTurb.r()[i] = vector(
							std::stod(aTurbConditions[3 * j]),
							std::stod(aTurbConditions[1 + 3 * j]),
							std::stod(aTurbConditions[2 + 3 * j]));
					break;
				}
			}

			if (!zoneFounded)
				throw exception(
						"Cell " + std::to_string(i)
								+ " of a dropped out of all zones.",
						errors::initialisationError);
		}

		/*Set initial conditions for b.*/
		for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
		{
			const vector & radiusOfCell { meshReference.cells()[i].rC() };
			bool zoneFounded { false };

			for (std::size_t j = 0; j < numberOfZones.size(); ++j)
			{
				zoneFounded = numberOfZones[j]->isPointInZone(radiusOfCell);

				if (zoneFounded)
				{
					phase->bTurb.r()[i] = std::stod(bTurbConditions[j]);
					break;
				}
			}

			if (!zoneFounded)
				throw exception(
						"Cell " + std::to_string(i)
								+ " of b dropped out of all zones.",
						errors::initialisationError);
		}
	}
	else
	{
		const auto noutputStr = std::to_string(readDataPoint.first);
		const std::size_t length_of_noutput = noutputStr.size();

		if (lengthOfNumber < length_of_noutput)
			throw exception(
					"Output number length larger than specified number of digits.",
					errors::tooBigOutputNumberError);

		std::string bufOutputN(lengthOfNumber - length_of_noutput, '0');

		bufOutputN.append(noutputStr);

		std::string fieldDataDirectoryName { "./fieldsOutput/output_" };
		fieldDataDirectoryName.append(bufOutputN);

		if (!std::filesystem::exists(fieldDataDirectoryName))
			throw exception(
					std::string("Directory ") + fieldDataDirectoryName
							+ std::string(" doesn't exits."),
					errors::systemError);

		fieldDataDirectoryName.append("/");

		const auto fieldDataFileName_a = fieldDataDirectoryName
				+ std::string("a.dat");
		const auto fieldDataFileName_b = fieldDataDirectoryName
				+ std::string("b.dat");
		const auto fieldDataFileName_epsilon = fieldDataDirectoryName
				+ std::string("epsilon.dat");
		const auto fieldDataFileName_k = fieldDataDirectoryName
				+ std::string("k.dat");
		const auto fieldDataFileName_pressure = fieldDataDirectoryName
				+ std::string("pressure.dat");
		const auto fieldDataFileName_velocity = fieldDataDirectoryName
				+ std::string("velocity.dat");
		std::vector<std::string> fieldDataFileName_concentration(
				phase->concentration.v.size() - 1);
		for (std::size_t k = 0; k < fieldDataFileName_concentration.size(); ++k)
			fieldDataFileName_concentration[k] = fieldDataDirectoryName
					+ std::string("sub_") + std::string(std::to_string(k + 1))
					+ std::string(".dat");

#ifdef MPI_VERSION
		const auto localCellInterval = std::make_pair(
				parallelism.mpi_rank * meshReference.cellsSize(),
				(parallelism.mpi_rank + 1) * meshReference.cellsSize());

		{
			std::ifstream input_aVector { fieldDataFileName_a };
			if (!input_aVector.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_a)
								+ std::string("."));
			input_aVector.precision(ioPrecision);

			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				std::string vectorData1, vectorData2, vectorData3;

				input_aVector >> vectorData1 >> vectorData2 >> vectorData3;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->aTurb.r()[i_local] = vector(std::stod(vectorData1),
							std::stod(vectorData2), std::stod(vectorData3));
				}
			}
			input_aVector.close();
		}

		{
			std::ifstream input_bScalar { fieldDataFileName_b };
			if (!input_bScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_b)
								+ std::string("."));
			input_bScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar bData;
				input_bScalar >> bData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->bTurb.r()[i_local] = bData;
				}
			}
			input_bScalar.close();
		}

		{
			std::ifstream input_epsilonScalar { fieldDataFileName_epsilon };
			if (!input_epsilonScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_epsilon)
								+ std::string("."));
			input_epsilonScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar epsilonData;
				input_epsilonScalar >> epsilonData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->epsTurb.r()[i_local] = epsilonData;
				}
			}
			input_epsilonScalar.close();
		}

		{
			std::ifstream input_kScalar { fieldDataFileName_k };
			if (!input_kScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_k)
								+ std::string("."));
			input_kScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar kData;
				input_kScalar >> kData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->kTurb.r()[i_local] = kData;
				}
			}
			input_kScalar.close();
		}

		{
			std::ifstream input_pressureScalar { fieldDataFileName_pressure };
			if (!input_pressureScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_pressure)
								+ std::string("."));
			input_pressureScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar pressureData;
				input_pressureScalar >> pressureData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->pressure.r()[i_local] = pressureData;
				}
			}
			input_pressureScalar.close();
		}

		{
			std::ifstream input_velocityVector { fieldDataFileName_velocity };
			if (!input_velocityVector.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_velocity)
								+ std::string("."));
			input_velocityVector.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				std::string vectorData1, vectorData2, vectorData3;

				input_velocityVector >> vectorData1 >> vectorData2
						>> vectorData3;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->velocity.r()[i_local] = vector(
							std::stod(vectorData1), std::stod(vectorData2),
							std::stod(vectorData3));
				}
			}
			input_velocityVector.close();
		}

		for (std::size_t k = 0; k < fieldDataFileName_concentration.size(); ++k)
		{
			std::ifstream input_concentrationScalar {
					fieldDataFileName_concentration[k] };
			if (!input_concentrationScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(
										fieldDataFileName_concentration[k])
								+ std::string("."));
			input_concentrationScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar concentrationData;
				input_concentrationScalar >> concentrationData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->concentration.v[k + 1].r()[i_local] =
							concentrationData;
				}
			}
			input_concentrationScalar.close();
		}
#else
		{
			std::ifstream input_aVector { fieldDataFileName_a };
			if (!input_aVector.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_a)
								+ std::string("."));
			input_aVector.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				std::string vectorData1, vectorData2, vectorData3;

				input_aVector >> vectorData1 >> vectorData2 >> vectorData3;

				phase->aTurb.r()[i] = vector(std::stod(vectorData1),
						std::stod(vectorData2), std::stod(vectorData3));
			}
			input_aVector.close();
		}

		{
			std::ifstream input_bScalar { fieldDataFileName_b };
			if (!input_bScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_b)
								+ std::string("."));
			input_bScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar bData;
				input_bScalar >> bData;

				phase->bTurb.r()[i] = bData;
			}
			input_bScalar.close();
		}

		{
			std::ifstream input_epsilonScalar { fieldDataFileName_epsilon };
			if (!input_epsilonScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_epsilon)
								+ std::string("."));
			input_epsilonScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar epsilonData;
				input_epsilonScalar >> epsilonData;

				phase->epsTurb.r()[i] = epsilonData;
			}
			input_epsilonScalar.close();
		}

		{
			std::ifstream input_kScalar { fieldDataFileName_k };
			if (!input_kScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_k)
								+ std::string("."));
			input_kScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar kData;
				input_kScalar >> kData;

				phase->kTurb.r()[i] = kData;
			}
			input_kScalar.close();
		}

		{
			std::ifstream input_pressureScalar { fieldDataFileName_pressure };
			if (!input_pressureScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_pressure)
								+ std::string("."));
			input_pressureScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar pressureData;
				input_pressureScalar >> pressureData;

				phase->pressure.r()[i] = pressureData;
			}
			input_pressureScalar.close();
		}

		{
			std::ifstream input_velocityVector { fieldDataFileName_velocity };
			if (!input_velocityVector.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_velocity)
								+ std::string("."));
			input_velocityVector.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				std::string vectorData1, vectorData2, vectorData3;

				input_velocityVector >> vectorData1 >> vectorData2
						>> vectorData3;

				phase->velocity.r()[i] = vector(std::stod(vectorData1),
						std::stod(vectorData2), std::stod(vectorData3));
			}
			input_velocityVector.close();
		}

		for (std::size_t k = 0; k < fieldDataFileName_concentration.size(); ++k)
		{
			std::ifstream input_concentrationScalar {
					fieldDataFileName_concentration[k] };
			if (!input_concentrationScalar.is_open())
				throw std::ifstream::failure(
						std::string("Couldn't open output file for field data ")
								+ std::string(
										fieldDataFileName_concentration[k])
								+ std::string("."));
			input_concentrationScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar concentrationData;
				input_concentrationScalar >> concentrationData;

				phase->concentration.v[k + 1].r()[i] = concentrationData;
			}
			input_concentrationScalar.close();
		}
#endif
	}

	/*Additional peak of k and epsilon*/
	{
		std::ifstream turbPeakConditionsFile { "./set/turbPeak.txt" };

		std::string profileType;

		scalar kMax, epsMax, xCenter, lBound, rBound;

		if (turbPeakConditionsFile.is_open())
			std::cout << "./set/turbPeak.txt is opened." << std::endl;
		else
			throw std::ifstream::failure("./set/turbPeak.txt not found.");

		turbPeakConditionsFile >> skipBuffer >> profileType >> skipBuffer
				>> kMax >> skipBuffer >> epsMax >> skipBuffer >> xCenter
				>> skipBuffer >> lBound >> skipBuffer >> rBound;

		turbPeakConditionsFile.close();

		skipBuffer.clear();

		if (profileType == "linear")
		{
			const auto leftX = xCenter - lBound;
			const auto rightX = xCenter + rBound;

			volumeField<scalar> kAdd { meshReference, 0 }, epsAdd {
					meshReference, 0 };

			const auto fL_k = kMax / lBound;
			const auto fR_k = kMax / rBound;
			const auto fL_eps = epsMax / lBound;
			const auto fR_eps = epsMax / rBound;

			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				const vector & cellRad = meshReference.cells()[i].rC();
				const scalar xCoord = std::get<0>(cellRad());
				if ((xCoord >= leftX) && (xCoord < xCenter))
				{
					const scalar relX = xCoord - leftX;
					kAdd.r()[i] = fL_k * relX;
					epsAdd.r()[i] = fL_eps * relX;
				}
				else if ((xCoord >= xCenter) && (xCoord <= rightX))
				{
					const scalar relX = xCoord - xCenter;
					kAdd.r()[i] = kMax - fR_k * relX;
					epsAdd.r()[i] = epsMax - fR_eps * relX;
				}
			}

			epsAdd.r() = std::pow(kAdd(), 3. / 2.) / (lBound + rBound);

			phase->kTurb.r() += kAdd();
			phase->epsTurb.r() += epsAdd();

			std::cout << "Added linear peak for k and epsilon." << std::endl;
		}
		else if (profileType != "no")
			throw exception(
					"Wrong parameter of peak, must be <<linear>> or <<no>>.",
					errors::initialisationError);
	}
	/**/

	/*Check for minimum values.*/
	std::replace_if(std::begin(phase->kTurb.r()), std::end(phase->kTurb.r()),
			[&phase](const scalar value) 
			{
				return value < phase->turbulenceSources->turbPar->mink();
			}, phase->turbulenceSources->turbPar->mink());

	std::replace_if(std::begin(phase->epsTurb.r()),
			std::end(phase->epsTurb.r()), [&phase](const scalar value) 
			{
				return value < phase->turbulenceSources->turbPar->mineps();
			}, phase->turbulenceSources->turbPar->mineps());

	/*Calculate derived fields.*/
	for (std::size_t k = 1; k < phase->concentration.v.size(); ++k)
		phase->concentration.v[0].r() += phase->concentration.v[k]();

	for (std::size_t k = 1; k < phase->density.size(); ++k)
		phase->density[k].r() = phase->concentration.v[k]()
				* phase->phaseThermodynamics->Mv()[k - 1];

	for (std::size_t k = 1; k < phase->density.size(); ++k)
		phase->density[0].r() += phase->density[k]();

	phase->momentum.r() = astProduct(phase->velocity, phase->density[0])();

	phase->internalEnergy.r() = phase->phaseThermodynamics->UvFromp(
			phase->concentration.p, phase->pressure());

	phase->temperature.r() = phase->phaseThermodynamics->TFromUv(
			phase->concentration.p, phase->internalEnergy());

	phase->rhokTurb.r() = phase->kTurb() * phase->density[0]();

	phase->rhoepsTurb.r() = phase->density[0]() * phase->epsTurb();

	phase->rhoaTurb.r() = astProduct(phase->aTurb, phase->density[0])();

	phase->rhobTurb.r() = phase->density[0]() * phase->bTurb();

	phase->calculateCoefficients(phase->kTurb, phase->epsTurb,
			*phase->turbulenceSources->turbPar);

	/*Set turbulent quantities to zero, if it is non-turbulent task.*/
	if (!phase->turbulenceSources->turbulence)
	{
		phase->kTurb.r() = scalar { 0 };
		phase->epsTurb.r() = 0;
		phase->rhokTurb.r() = scalar { 0 };
		phase->rhoepsTurb.r() = 0;
		phase->aTurb.r() = vector(0);
		phase->bTurb.r() = 0;
		phase->rhoaTurb.r() = vector(0);
		phase->rhobTurb.r() = 0;

		phase->tAssign(0.);
	}

	{
		const auto v2 = ampProduct(phase->velocity, phase->velocity);

		phase->totalEnergy.r() = phase->internalEnergy()
				+ phase->density[0]() * v2() * 0.5 + phase->rhokTurb();
	}

	phase->HelmholtzEnergy.r() = phase->phaseThermodynamics->Fv(
			phase->concentration.p, phase->temperature());

	phase->entropy.r() = phase->phaseThermodynamics->Sv(phase->concentration.p,
			phase->temperature());

	return std::make_tuple(std::move(phase), enthalpyFlowFlag,
			molMassDiffusionFlag);
}
