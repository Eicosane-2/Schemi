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
#include "fabricFunctions.hpp"
#include "mixtureIdeal.hpp"
#include "mixtureRedlichKwong.hpp"
#include "mixtureVanDerWaals.hpp"
#include "fieldProducts.hpp"

std::tuple<std::unique_ptr<schemi::homogeneousPhase<schemi::cubicCell>>,
		schemi::enthalpyFlowEnum, bool> schemi::phaseInitialization(
		std::size_t numberOfComponents, std::size_t numberOfZones,
		const mesh & meshReference,
		const std::vector<boundaryConditionType> & commonConditions,
		const MPIHandler & parallelism, const std::string & turbulenceONString,
		const std::string & sourceTypeString,
		const std::string & universalGasConstant,
		const std::string & equationOfState, const std::size_t readDataPoint)
{
	std::unique_ptr<homogeneousPhase<cubicCell>> phase { nullptr };
	enthalpyFlowEnum enthalpyFlowFlag;
	bool molMassDiffusionFlag;

	std::string skipBuffer;

	std::unique_ptr<abstractTurbulenceGen> turbulenceSources(
			createTurbulenceModel(turbulenceONString, sourceTypeString));

	transportCoefficients<cubicCell> tCoeffsPhase { meshReference };

	bunchOfFields<cubicCell> cellFields { meshReference, numberOfComponents };
	std::unique_ptr<abstractMixtureThermodynamics> mixture;

	/*Set boundary condition for transport coefficients.*/
	auto bndConCom = commonConditions;

	tCoeffsPhase.tD = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });
	tCoeffsPhase.tKappa = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });
	tCoeffsPhase.tLambda = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });
	tCoeffsPhase.k_D = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });
	tCoeffsPhase.eps_D = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });
	tCoeffsPhase.a_D = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });
	tCoeffsPhase.b_D = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });

	std::replace(bndConCom.begin(), bndConCom.end(),
			boundaryConditionType::calculated,
			boundaryConditionType::calculatedTurbulentViscosity);
	tCoeffsPhase.tNu = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });

	std::ifstream transportCoefficientsFile { "./set/transportCoefficients.txt",
			std::ios::in };
	if (transportCoefficientsFile.is_open())
		std::cout << "./set/transportCoefficients.txt is opened." << std::endl;
	else
		throw exception("./set/transportCoefficients.txt not found.",
				errorsEnum::initializationError);

	bndConCom = commonConditions;
	std::replace(bndConCom.begin(), bndConCom.end(),
			boundaryConditionType::calculated,
			boundaryConditionType::freeBoundary);
	tCoeffsPhase.pNu = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });
	tCoeffsPhase.pD = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });
	tCoeffsPhase.pKappa = volumeField<scalar>(meshReference, 0,
			subPatchData<scalar> { bndConCom[0] }, subPatchData<scalar> {
					bndConCom[1] }, subPatchData<scalar> { bndConCom[2] },
			subPatchData<scalar> { bndConCom[3] }, subPatchData<scalar> {
					bndConCom[4] }, subPatchData<scalar> { bndConCom[5] });
	scalar physNu, physD, physKappa;
	std::string implicitEnthalpyFlowFlagString, MolMassDiffusionFlagString;
	transportCoefficientsFile >> skipBuffer >> physNu >> skipBuffer >> physD
			>> skipBuffer >> physKappa >> skipBuffer
			>> implicitEnthalpyFlowFlagString >> skipBuffer
			>> MolMassDiffusionFlagString;
	tCoeffsPhase.pNu.ref_r() = physNu;
	tCoeffsPhase.pD.ref_r() = physD;
	tCoeffsPhase.pKappa.ref_r() = physKappa;
	transportCoefficientsFile.close();

	if (implicitEnthalpyFlowFlagString == "implicit")
		enthalpyFlowFlag = enthalpyFlowEnum::implicitSolve;
	else if (implicitEnthalpyFlowFlagString == "explicit")
		enthalpyFlowFlag = enthalpyFlowEnum::explicitSolve;
	else if (implicitEnthalpyFlowFlagString == "no")
		enthalpyFlowFlag = enthalpyFlowEnum::noSolve;
	else
		throw exception("Unknown flag for enthalpy flow calculation.",
				errorsEnum::initializationError);

	if (MolMassDiffusionFlagString == "on")
		molMassDiffusionFlag = true;
	else if (MolMassDiffusionFlagString == "off")
		molMassDiffusionFlag = false;
	else
		throw exception("Unknown flag for molar mass diffusion correction.",
				errorsEnum::initializationError);

	/*Create mixture class and boundary conditions for basic quantities: concentrations, velocity, pressure, k, epsilon, a and b.*/
	std::vector<std::array<std::vector<subPatchData<scalar>>, 6>> boundaryConditionsMatrix(
			numberOfComponents);
	for (auto & arr : boundaryConditionsMatrix)
		arr.fill(std::vector<subPatchData<scalar>>(0));

	std::vector<std::vector<std::string>> matrixOfSubstancesConditions(
			numberOfComponents, std::vector<std::string>(4 + numberOfZones));

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		std::string substanceName { "./set/sub_" }, bufComponent(
				std::to_string(k + 1));

		substanceName.append(bufComponent);
		substanceName.append(".txt");
		std::ifstream substanceConditionsFile { substanceName, std::ios::in };

		if (substanceConditionsFile.is_open())
			std::cout << substanceName << " is opened." << std::endl;
		else
			throw exception(substanceName + " not found.",
					errorsEnum::initializationError);

		substanceConditionsFile >> skipBuffer
				>> matrixOfSubstancesConditions[k][0] /*M*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][1] /*Cv*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][2] /*Tcrit*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][3]; /*Pcrit*/

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

				if (boundTp == boundaryConditionType::fixedValue)
				{
					scalar fv;
					substanceConditionsFile >> fv;

					boundaryConditionsMatrix[k][boundaryI].push_back( { boundTp,
							fv });
				}
				else
					boundaryConditionsMatrix[k][boundaryI].push_back(
							{ boundTp });
			}
			else
			{
				std::size_t subPatchNumber;

				substanceConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errorsEnum::initializationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					substanceConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValue)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						substanceConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						substanceConditionsFile >> fv;

						boundaryConditionsMatrix[k][boundaryI].push_back( {
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 }, fv });
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						substanceConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsMatrix[k][boundaryI].push_back( {
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 } });
					}
				}
			}
		}

		substanceConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errorsEnum::initializationError);
		for (std::size_t j = 0; j < numberOfZones; ++j)
			substanceConditionsFile >> matrixOfSubstancesConditions[k][4 + j]; /*Values in zones*/
		substanceConditionsFile.close();
	}

	/*Convert strings to values.*/
	std::array<std::valarray<scalar>, 4> thermodynamicalProperties;
	thermodynamicalProperties.fill(std::valarray<scalar>(numberOfComponents));

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		thermodynamicalProperties[0][k] = std::stod(
				matrixOfSubstancesConditions[k][0]); /*M*/
		thermodynamicalProperties[1][k] = std::stod(
				matrixOfSubstancesConditions[k][1]); /*Cv*/
		thermodynamicalProperties[2][k] = std::stod(
				matrixOfSubstancesConditions[k][2]); /*Tcrit*/
		thermodynamicalProperties[3][k] = std::stod(
				matrixOfSubstancesConditions[k][3]); /*Pcrit*/
	}

	for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
	{
		parallelism.correctBoundaryConditions(boundaryConditionsMatrix[k - 1]);

		cellFields.concentration.v[k] = volumeField<scalar>(meshReference, 0,
				boundaryConditionsMatrix[k - 1][0],
				boundaryConditionsMatrix[k - 1][1],
				boundaryConditionsMatrix[k - 1][2],
				boundaryConditionsMatrix[k - 1][3],
				boundaryConditionsMatrix[k - 1][4],
				boundaryConditionsMatrix[k - 1][5]);

		auto bndVConRho_k = boundaryConditionsMatrix[k - 1];

		for (auto & b_s : bndVConRho_k)
			for (auto & b_i : b_s)
				b_i.fixVal *= thermodynamicalProperties[0][k - 1];

		cellFields.density[k] = volumeField<scalar>(meshReference, 0,
				bndVConRho_k[0], bndVConRho_k[1], bndVConRho_k[2],
				bndVConRho_k[3], bndVConRho_k[4], bndVConRho_k[5]);
	}
	{
		cellFields.concentration.v[0] = volumeField<scalar>(meshReference, 0);

		auto bndConRho = commonConditions;

		std::replace(bndConRho.begin(), bndConRho.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedDensity);

		cellFields.density[0] = volumeField<scalar>(meshReference, 0,
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
				errorsEnum::initializationError);

	{
		gasModelEnum gasModelFlag;
		if (equationOfState == "vanDerWaals")
			gasModelFlag = gasModelEnum::vanDerWaals;
		else if (equationOfState == "ideal")
			gasModelFlag = gasModelEnum::ideal;
		else if (equationOfState == "RedlichKwong")
			gasModelFlag = gasModelEnum::RedlichKwong;
		else
			throw exception("Unknown type of equation of state.",
					errorsEnum::initializationError);

		switch (gasModelFlag)
		{
		case gasModelEnum::vanDerWaals:
			mixture = std::make_unique<mixtureVanDerWaals>(R, hPlanck,
					thermodynamicalProperties[0], thermodynamicalProperties[1],
					thermodynamicalProperties[2], thermodynamicalProperties[3]);
			break;
		case gasModelEnum::RedlichKwong:
			mixture = std::make_unique<mixtureRedlichKwong>(R, hPlanck,
					thermodynamicalProperties[0], thermodynamicalProperties[1],
					thermodynamicalProperties[2], thermodynamicalProperties[3]);
			break;
		case gasModelEnum::ideal:
		default:
			mixture = std::make_unique<mixtureIdeal>(R, hPlanck,
					thermodynamicalProperties[0], thermodynamicalProperties[1],
					thermodynamicalProperties[2], thermodynamicalProperties[3]);
			break;
		}
	}

	phase = std::make_unique<homogeneousPhase<cubicCell>>(cellFields,
			tCoeffsPhase, mixture, turbulenceSources);

	std::array<std::vector<subPatchData<vector>>, 6> boundaryConditionsVelocity;
	boundaryConditionsVelocity.fill(std::vector<subPatchData<vector>>(0));
	std::vector<std::string> velocityConditions(3 * numberOfZones);
	{
		std::ifstream velocityConditionsFile { "./set/velocity.txt",
				std::ios::in };

		if (velocityConditionsFile.is_open())
			std::cout << "./set/velocity.txt is opened." << std::endl;
		else
			throw exception("./set/velocity.txt not found.",
					errorsEnum::initializationError);

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

				if (boundTp == boundaryConditionType::fixedValue)
				{
					scalar fv1, fv2, fv3;
					velocityConditionsFile >> fv1 >> fv2 >> fv3;

					boundaryConditionsVelocity[boundaryI].push_back( { boundTp,
							vector { fv1, fv2, fv3 } });
				}
				else
					boundaryConditionsVelocity[boundaryI].push_back(
							{ boundTp });
			}
			else
			{
				std::size_t subPatchNumber;

				velocityConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errorsEnum::initializationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					velocityConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValue)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						velocityConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv1, fv2, fv3;
						velocityConditionsFile >> fv1 >> fv2 >> fv3;

						boundaryConditionsVelocity[boundaryI].push_back( {
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 }, vector { fv1, fv2, fv3 } });
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						velocityConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsVelocity[boundaryI].push_back( {
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 } });
					}
				}
			}
		}

		velocityConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errorsEnum::initializationError);
		for (std::size_t i = 0; i < numberOfZones; ++i)
			velocityConditionsFile >> velocityConditions[3 * i]
					>> velocityConditions[1 + 3 * i]
					>> velocityConditions[2 + 3 * i];

		velocityConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsVelocity);

		phase->velocity = volumeField<vector>(meshReference, vector(0),
				boundaryConditionsVelocity[0], boundaryConditionsVelocity[1],
				boundaryConditionsVelocity[2], boundaryConditionsVelocity[3],
				boundaryConditionsVelocity[4], boundaryConditionsVelocity[5]);
	}

	std::array<std::vector<subPatchData<scalar>>, 6> boundaryConditionsPressure;
	boundaryConditionsPressure.fill(std::vector<subPatchData<scalar>>(0));
	std::vector<std::string> pressureConditions(numberOfZones);
	{
		std::ifstream pressureConditionsFile { "./set/pressure.txt",
				std::ios::in };

		if (pressureConditionsFile.is_open())
			std::cout << "./set/pressure.txt is opened." << std::endl;
		else
			throw exception("./set/pressure.txt not found.",
					errorsEnum::initializationError);

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

				if (boundTp == boundaryConditionType::fixedValue)
				{
					scalar fv;
					pressureConditionsFile >> fv;

					boundaryConditionsPressure[boundaryI].push_back( { boundTp,
							fv });
				}
				else
					boundaryConditionsPressure[boundaryI].push_back(
							{ boundTp });
			}
			else
			{
				std::size_t subPatchNumber;

				pressureConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errorsEnum::initializationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					pressureConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValue)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						pressureConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						pressureConditionsFile >> fv;

						boundaryConditionsPressure[boundaryI].push_back( {
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 }, fv });
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						pressureConditionsFile >> skipBuffer >> vb1 >> vb2
								>> vb3 >> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsPressure[boundaryI].push_back( {
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 } });
					}
				}
			}
		}

		pressureConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errorsEnum::initializationError);
		for (std::size_t i = 0; i < numberOfZones; ++i)
			pressureConditionsFile >> pressureConditions[i];

		pressureConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsPressure);

		phase->pressure = volumeField<scalar>(meshReference, 0,
				boundaryConditionsPressure[0], boundaryConditionsPressure[1],
				boundaryConditionsPressure[2], boundaryConditionsPressure[3],
				boundaryConditionsPressure[4], boundaryConditionsPressure[5]);
	}

	std::array<std::vector<subPatchData<scalar>>, 6> boundaryConditionskTurb;
	boundaryConditionskTurb.fill(std::vector<subPatchData<scalar>>(0));
	std::vector<std::string> kTurbConditions(numberOfZones);
	{
		std::ifstream kTurbConditionsFile { "./set/k.txt", std::ios::in };

		if (kTurbConditionsFile.is_open())
			std::cout << "./set/k.txt is opened." << std::endl;
		else
			throw exception("./set/k.txt is opened.",
					errorsEnum::initializationError);

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

				if (boundTp == boundaryConditionType::fixedValue)
				{
					scalar fv;
					kTurbConditionsFile >> fv;

					boundaryConditionskTurb[boundaryI].push_back(
							{ boundTp, fv });
				}
				else
					boundaryConditionskTurb[boundaryI].push_back( { boundTp });
			}
			else
			{
				std::size_t subPatchNumber;

				kTurbConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errorsEnum::initializationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					kTurbConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValue)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						kTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						kTurbConditionsFile >> fv;

						boundaryConditionskTurb[boundaryI].push_back( { boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 }, fv });
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						kTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionskTurb[boundaryI].push_back( { boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 } });
					}
				}
			}
		}

		kTurbConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errorsEnum::initializationError);
		for (std::size_t i = 0; i < numberOfZones; ++i)
			kTurbConditionsFile >> kTurbConditions[i];

		kTurbConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionskTurb);

		phase->kTurb = volumeField<scalar>(meshReference, scalar { 0 },
				boundaryConditionskTurb[0], boundaryConditionskTurb[1],
				boundaryConditionskTurb[2], boundaryConditionskTurb[3],
				boundaryConditionskTurb[4], boundaryConditionskTurb[5]);
	}

	std::array<std::vector<subPatchData<scalar>>, 6> boundaryConditionsepsTurb;
	boundaryConditionsepsTurb.fill(std::vector<subPatchData<scalar>>(0));
	std::vector<std::string> epsTurbConditions(numberOfZones);
	{
		std::ifstream epsTurbConditionsFile { "./set/epsilon.txt", std::ios::in };

		if (epsTurbConditionsFile.is_open())
			std::cout << "./set/epsilon.txt is opened." << std::endl;
		else
			throw exception("./set/epsilon.txt not found.",
					errorsEnum::initializationError);

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

				if (boundTp == boundaryConditionType::fixedValue)
				{
					scalar fv;
					epsTurbConditionsFile >> fv;

					boundaryConditionsepsTurb[boundaryI].push_back( { boundTp,
							fv });
				}
				else
					boundaryConditionsepsTurb[boundaryI].push_back(
							{ boundTp });
			}
			else
			{
				std::size_t subPatchNumber;

				epsTurbConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errorsEnum::initializationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					epsTurbConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValue)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						epsTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						epsTurbConditionsFile >> fv;

						boundaryConditionsepsTurb[boundaryI].push_back( {
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 }, fv });
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						epsTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsepsTurb[boundaryI].push_back( {
								boundTp, vector { vb1, vb2, vb3 }, vector { ve1,
										ve2, ve3 } });
					}
				}
			}
		}

		epsTurbConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errorsEnum::initializationError);
		for (std::size_t i = 0; i < numberOfZones; ++i)
			epsTurbConditionsFile >> epsTurbConditions[i];

		epsTurbConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsepsTurb);

		phase->epsTurb = volumeField<scalar>(meshReference, 0,
				boundaryConditionsepsTurb[0], boundaryConditionsepsTurb[1],
				boundaryConditionsepsTurb[2], boundaryConditionsepsTurb[3],
				boundaryConditionsepsTurb[4], boundaryConditionsepsTurb[5]);
	}

	std::array<std::vector<subPatchData<vector>>, 6> boundaryConditionsaTurb;
	boundaryConditionsaTurb.fill(std::vector<subPatchData<vector>>(0));
	std::vector<std::string> aTurbConditions(3 * numberOfZones);
	{
		std::ifstream aTurbConditionsFile { "./set/a.txt", std::ios::in };

		if (aTurbConditionsFile.is_open())
			std::cout << "./set/a.txt is opened." << std::endl;
		else
			throw exception("./set/a.txt not found.",
					errorsEnum::initializationError);

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

				if (boundTp == boundaryConditionType::fixedValue)
				{
					scalar fv1, fv2, fv3;
					aTurbConditionsFile >> fv1 >> fv2 >> fv3;

					boundaryConditionsaTurb[boundaryI].push_back( { boundTp,
							vector { fv1, fv2, fv3 } });
				}
				else
					boundaryConditionsaTurb[boundaryI].push_back( { boundTp });
			}
			else
			{
				std::size_t subPatchNumber;

				aTurbConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errorsEnum::initializationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					aTurbConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValue)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						aTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv1, fv2, fv3;
						aTurbConditionsFile >> fv1 >> fv2 >> fv3;

						boundaryConditionsaTurb[boundaryI].push_back( { boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 },
								vector { fv1, fv2, fv3 } });
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						aTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsaTurb[boundaryI].push_back( { boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 } });
					}
				}
			}
		}

		aTurbConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errorsEnum::initializationError);
		for (std::size_t i = 0; i < numberOfZones; ++i)
			aTurbConditionsFile >> aTurbConditions[3 * i]
					>> aTurbConditions[1 + 3 * i] >> aTurbConditions[2 + 3 * i];

		aTurbConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsaTurb);

		phase->aTurb = volumeField<vector>(meshReference, vector(0, 0, 0),
				boundaryConditionsaTurb[0], boundaryConditionsaTurb[1],
				boundaryConditionsaTurb[2], boundaryConditionsaTurb[3],
				boundaryConditionsaTurb[4], boundaryConditionsaTurb[5]);
	}

	std::array<std::vector<subPatchData<scalar>>, 6> boundaryConditionsbTurb;
	boundaryConditionsbTurb.fill(std::vector<subPatchData<scalar>>(0));
	std::vector<std::string> bTurbConditions(numberOfZones);
	{
		std::ifstream bTurbConditionsFile { "./set/b.txt", std::ios::in };

		if (bTurbConditionsFile.is_open())
			std::cout << "./set/b.txt is opened." << std::endl;
		else
			throw exception("./set/b.txt not found.",
					errorsEnum::initializationError);

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

				if (boundTp == boundaryConditionType::fixedValue)
				{
					scalar fv;
					bTurbConditionsFile >> fv;

					boundaryConditionsbTurb[boundaryI].push_back(
							{ boundTp, fv });
				}
				else
					boundaryConditionsbTurb[boundaryI].push_back( { boundTp });
			}
			else
			{
				std::size_t subPatchNumber;

				bTurbConditionsFile >> subPatchNumber;

				if (subPatchNumber == 1)
					throw exception("SubPatches must be more than one.",
							errorsEnum::initializationError);

				for (std::size_t sp = 0; sp < subPatchNumber; ++sp)
				{
					std::string subPatchBoundary;
					bTurbConditionsFile >> subPatchBoundary;

					const auto boundTp = boundaryConditionFromString(
							subPatchBoundary);

					if (boundTp == boundaryConditionType::fixedValue)
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						bTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						scalar fv;
						bTurbConditionsFile >> fv;

						boundaryConditionsbTurb[boundaryI].push_back( { boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 }, fv });
					}
					else
					{
						scalar vb1, vb2, vb3, ve1, ve2, ve3;

						bTurbConditionsFile >> skipBuffer >> vb1 >> vb2 >> vb3
								>> skipBuffer >> ve1 >> ve2 >> ve3
								>> skipBuffer;

						boundaryConditionsbTurb[boundaryI].push_back( { boundTp,
								vector { vb1, vb2, vb3 },
								vector { ve1, ve2, ve3 } });
					}
				}
			}
		}

		bTurbConditionsFile >> skipBuffer;
		if (skipBuffer != "#Values_in_zones")
			throw exception("Wrong position in file.",
					errorsEnum::initializationError);
		for (std::size_t i = 0; i < numberOfZones; ++i)
			bTurbConditionsFile >> bTurbConditions[i];

		bTurbConditionsFile.close();

		parallelism.correctBoundaryConditions(boundaryConditionsbTurb);

		phase->bTurb = volumeField<scalar>(meshReference,
				phase->turbulenceSources->turbPar->minb_value,
				boundaryConditionsbTurb[0], boundaryConditionsbTurb[1],
				boundaryConditionsbTurb[2], boundaryConditionsbTurb[3],
				boundaryConditionsbTurb[4], boundaryConditionsbTurb[5]);
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

	/*Reading vectors for zones determination.*/
	std::valarray<vector> vectorsOfZones(numberOfZones * 2);
	{
		std::ifstream coordinatesOfZonesFile("./set/coordinatesOfZones.txt",
				std::ios::in);

		if (coordinatesOfZonesFile.is_open())
			std::cout << "./set/coordinatesOfZones.txt is opened." << std::endl;
		else
			throw exception("./set/coordinatesOfZones.txt not found.",
					errorsEnum::initializationError);

		for (std::size_t i = 0; i < numberOfZones; ++i)
		{
			coordinatesOfZonesFile >> skipBuffer
					>> vectorsOfZones[2 * i].v_r()[0]
					>> vectorsOfZones[2 * i].v_r()[1]
					>> vectorsOfZones[2 * i].v_r()[2];
			coordinatesOfZonesFile >> skipBuffer
					>> vectorsOfZones[2 * i + 1].v_r()[0]
					>> vectorsOfZones[2 * i + 1].v_r()[1]
					>> vectorsOfZones[2 * i + 1].v_r()[2];
		}

		coordinatesOfZonesFile.close();
	}

	/*Set initial conditions for concentrations.*/
	for (std::size_t k = 1; k < phase->concentration.v.size(); ++k)
		for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
		{
			const vector & radiusOfCell { meshReference.cells()[i].rC() };
			bool zoneFounded { false };

			for (std::size_t j = 0; j < numberOfZones; ++j)
			{
				if ((radiusOfCell.v()[0] > vectorsOfZones[2 * j].v()[0])
						&& (radiusOfCell.v()[0]
								< vectorsOfZones[2 * j + 1].v()[0])
						&& (radiusOfCell.v()[1] > vectorsOfZones[2 * j].v()[1])
						&& (radiusOfCell.v()[1]
								< vectorsOfZones[2 * j + 1].v()[1])
						&& (radiusOfCell.v()[2] > vectorsOfZones[2 * j].v()[2])
						&& (radiusOfCell.v()[2]
								< vectorsOfZones[2 * j + 1].v()[2]))
				{
					zoneFounded = true;
					phase->concentration.v[k].ref_r()[i] = std::stod(
							matrixOfSubstancesConditions[k - 1][4 + j]);
					break;
				}
			}

			if (!zoneFounded)
				throw exception(
						"Cell " + std::to_string(i) + " of " + std::to_string(k)
								+ " component concentration dropped out of all zones.",
						errorsEnum::initializationError);
		}

	/*Set initial conditions for velocity.*/
	for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
	{
		const vector & radiusOfCell { meshReference.cells()[i].rC() };
		bool zoneFounded { false };

		for (std::size_t j = 0; j < numberOfZones; ++j)
		{
			if ((radiusOfCell.v()[0] > vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] > vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] > vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->velocity.ref_r()[i] = vector(
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
					errorsEnum::initializationError);
	}

	/*Set initial conditions for pressure.*/
	for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
	{
		const vector & radiusOfCell { meshReference.cells()[i].rC() };
		bool zoneFounded { false };

		for (std::size_t j = 0; j < numberOfZones; ++j)
		{
			if ((radiusOfCell.v()[0] > vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] > vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] > vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->pressure.ref_r()[i] = std::stod(pressureConditions[j]);
				break;
			}
		}
		if (!zoneFounded)
			throw exception(
					"Cell " + std::to_string(i)
							+ " of pressure dropped out of all zones.",
					errorsEnum::initializationError);
	}

	/*Set initial conditions for k.*/
	for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
	{
		const vector & radiusOfCell { meshReference.cells()[i].rC() };
		bool zoneFounded { false };

		for (std::size_t j = 0; j < numberOfZones; ++j)
		{
			if ((radiusOfCell.v()[0] > vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] > vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] > vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->kTurb.ref_r()[i] = std::stod(kTurbConditions[j]);
				break;
			}
		}
		if (!zoneFounded)
			throw exception(
					"Cell " + std::to_string(i)
							+ " of k dropped out of all zones.",
					errorsEnum::initializationError);
	}

	/*Set initial conditions for epsilon.*/
	for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
	{
		const vector & radiusOfCell { meshReference.cells()[i].rC() };
		bool zoneFounded { false };

		for (std::size_t j = 0; j < numberOfZones; ++j)
		{
			if ((radiusOfCell.v()[0] > vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] > vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] > vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->epsTurb.ref_r()[i] = std::stod(epsTurbConditions[j]);
				break;
			}
		}
		if (!zoneFounded)
			throw exception(
					"Cell " + std::to_string(i)
							+ " of epsilon dropped out of all zones.",
					errorsEnum::initializationError);
	}

	/*Set initial conditions for a.*/
	for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
	{
		const vector & radiusOfCell { meshReference.cells()[i].rC() };
		bool zoneFounded { false };

		for (std::size_t j = 0; j < numberOfZones; ++j)
		{
			if ((radiusOfCell.v()[0] > vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] > vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] > vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->aTurb.ref_r()[i] = vector(
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
					errorsEnum::initializationError);
	}

	/*Set initial conditions for b.*/
	for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
	{
		const vector & radiusOfCell { meshReference.cells()[i].rC() };
		bool zoneFounded { false };

		for (std::size_t j = 0; j < numberOfZones; ++j)
		{
			if ((radiusOfCell.v()[0] > vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] > vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] > vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->bTurb.ref_r()[i] = std::stod(bTurbConditions[j]);
				break;
			}
		}
		if (!zoneFounded)
			throw exception(
					"Cell " + std::to_string(i)
							+ " of b dropped out of all zones.",
					errorsEnum::initializationError);
	}

	if (readDataPoint)
	{
		const auto noutputStr = std::to_string(readDataPoint);
		const std::size_t length_of_noutput = noutputStr.size();

		if (lengthOfNumber < length_of_noutput)
			throw exception(
					"Output number length larger than specified number of digits.",
					errorsEnum::tooBigOutputNumberError);

		std::string bufOutputN(lengthOfNumber - length_of_noutput, '0');

		bufOutputN.append(noutputStr);

		std::string fieldDataDirectoryName { "./fieldsOutput/output_" };
		fieldDataDirectoryName.append(bufOutputN);

		if (!std::filesystem::exists(fieldDataDirectoryName))
			throw exception(
					std::string("Directory ") + fieldDataDirectoryName
							+ std::string(" doesn't exits."),
					errorsEnum::systemError);

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
			std::ifstream input_aVector { fieldDataFileName_a, std::ios::in };
			if (!input_aVector.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_a)
								+ std::string("."), errorsEnum::systemError);
			input_aVector.precision(ioPrecision);

			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				std::string vectorData1, vectorData2, vectorData3;

				input_aVector >> vectorData1 >> vectorData2 >> vectorData3;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->aTurb.ref_r()[i_local] = vector(
							std::stod(vectorData1), std::stod(vectorData2),
							std::stod(vectorData3));
				}
			}
			input_aVector.close();
		}

		{
			std::ifstream input_bScalar { fieldDataFileName_b, std::ios::in };
			if (!input_bScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_b)
								+ std::string("."), errorsEnum::systemError);
			input_bScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar bData;
				input_bScalar >> bData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->bTurb.ref_r()[i_local] = bData;
				}
			}
			input_bScalar.close();
		}

		{
			std::ifstream input_epsilonScalar { fieldDataFileName_epsilon,
					std::ios::in };
			if (!input_epsilonScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_epsilon)
								+ std::string("."), errorsEnum::systemError);
			input_epsilonScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar epsilonData;
				input_epsilonScalar >> epsilonData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->epsTurb.ref_r()[i_local] = epsilonData;
				}
			}
			input_epsilonScalar.close();
		}

		{
			std::ifstream input_kScalar { fieldDataFileName_k, std::ios::in };
			if (!input_kScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_k)
								+ std::string("."), errorsEnum::systemError);
			input_kScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar kData;
				input_kScalar >> kData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->kTurb.ref_r()[i_local] = kData;
				}
			}
			input_kScalar.close();
		}

		{
			std::ifstream input_pressureScalar { fieldDataFileName_pressure,
					std::ios::in };
			if (!input_pressureScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_pressure)
								+ std::string("."), errorsEnum::systemError);
			input_pressureScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar pressureData;
				input_pressureScalar >> pressureData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->pressure.ref_r()[i_local] = pressureData;
				}
			}
			input_pressureScalar.close();
		}

		{
			std::ifstream input_velocityVector { fieldDataFileName_velocity,
					std::ios::in };
			if (!input_velocityVector.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_velocity)
								+ std::string("."), errorsEnum::systemError);
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

					phase->velocity.ref_r()[i_local] = vector(
							std::stod(vectorData1), std::stod(vectorData2),
							std::stod(vectorData3));
				}
			}
			input_velocityVector.close();
		}

		for (std::size_t k = 0; k < fieldDataFileName_concentration.size(); ++k)
		{
			std::ifstream input_concentrationScalar {
					fieldDataFileName_concentration[k], std::ios::in };
			if (!input_concentrationScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(
										fieldDataFileName_concentration[k])
								+ std::string("."), errorsEnum::systemError);
			input_concentrationScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < parallelism.totCellNum(); ++i)
			{
				scalar concentrationData;
				input_concentrationScalar >> concentrationData;

				if ((i >= localCellInterval.first)
						&& (i < localCellInterval.second))
				{
					const auto i_local = i - localCellInterval.first;

					phase->concentration.v[k + 1].ref_r()[i_local] =
							concentrationData;
				}
			}
			input_concentrationScalar.close();
		}
#else
		{
			std::ifstream input_aVector { fieldDataFileName_a, std::ios::in };
			if (!input_aVector.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_a)
								+ std::string("."), errorsEnum::systemError);
			input_aVector.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				std::string vectorData1, vectorData2, vectorData3;

				input_aVector >> vectorData1 >> vectorData2 >> vectorData3;

				phase->aTurb.ref_r()[i] = vector(std::stod(vectorData1),
						std::stod(vectorData2), std::stod(vectorData3));
			}
			input_aVector.close();
		}

		{
			std::ifstream input_bScalar { fieldDataFileName_b, std::ios::in };
			if (!input_bScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_b)
								+ std::string("."), errorsEnum::systemError);
			input_bScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar bData;
				input_bScalar >> bData;

				phase->bTurb.ref_r()[i] = bData;
			}
			input_bScalar.close();
		}

		{
			std::ifstream input_epsilonScalar { fieldDataFileName_epsilon,
					std::ios::in };
			if (!input_epsilonScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_epsilon)
								+ std::string("."), errorsEnum::systemError);
			input_epsilonScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar epsilonData;
				input_epsilonScalar >> epsilonData;

				phase->epsTurb.ref_r()[i] = epsilonData;
			}
			input_epsilonScalar.close();
		}

		{
			std::ifstream input_kScalar { fieldDataFileName_k, std::ios::in };
			if (!input_kScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_k)
								+ std::string("."), errorsEnum::systemError);
			input_kScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar kData;
				input_kScalar >> kData;

				phase->kTurb.ref_r()[i] = kData;
			}
			input_kScalar.close();
		}

		{
			std::ifstream input_pressureScalar { fieldDataFileName_pressure,
					std::ios::in };
			if (!input_pressureScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_pressure)
								+ std::string("."), errorsEnum::systemError);
			input_pressureScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar pressureData;
				input_pressureScalar >> pressureData;

				phase->pressure.ref_r()[i] = pressureData;
			}
			input_pressureScalar.close();
		}

		{
			std::ifstream input_velocityVector { fieldDataFileName_velocity,
					std::ios::in };
			if (!input_velocityVector.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(fieldDataFileName_velocity)
								+ std::string("."), errorsEnum::systemError);
			input_velocityVector.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				std::string vectorData1, vectorData2, vectorData3;

				input_velocityVector >> vectorData1 >> vectorData2
						>> vectorData3;

				phase->velocity.ref_r()[i] = vector(std::stod(vectorData1),
						std::stod(vectorData2), std::stod(vectorData3));
			}
			input_velocityVector.close();
		}

		for (std::size_t k = 0; k < fieldDataFileName_concentration.size(); ++k)
		{
			std::ifstream input_concentrationScalar {
					fieldDataFileName_concentration[k], std::ios::in };
			if (!input_concentrationScalar.is_open())
				throw exception(
						std::string("Couldn't open output file for field data ")
								+ std::string(
										fieldDataFileName_concentration[k])
								+ std::string("."), errorsEnum::systemError);
			input_concentrationScalar.precision(ioPrecision);
			for (std::size_t i = 0; i < meshReference.cellsSize(); ++i)
			{
				scalar concentrationData;
				input_concentrationScalar >> concentrationData;

				phase->concentration.v[k + 1].ref_r()[i] = concentrationData;
			}
			input_concentrationScalar.close();
		}
#endif
	}

	/*Additional peak of k and epsilon*/
	{
		std::ifstream turbPeakConditionsFile { "./set/turbPeak.txt",
				std::ios::in };

		std::string profileType;

		scalar kMax, epsMax, xCenter, lBound, rBound;

		if (turbPeakConditionsFile.is_open())
			std::cout << "./set/turbPeak.txt is opened." << std::endl;
		else
			throw exception("./set/turbPeak.txt not found.",
					errorsEnum::initializationError);

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
				const scalar xCoord = cellRad.v()[0];
				if ((xCoord >= leftX) && (xCoord < xCenter))
				{
					const scalar relX = xCoord - leftX;
					kAdd.ref_r()[i] = fL_k * relX;
					epsAdd.ref_r()[i] = fL_eps * relX;
				}
				else if ((xCoord >= xCenter) && (xCoord <= rightX))
				{
					const scalar relX = xCoord - xCenter;
					kAdd.ref_r()[i] = kMax - fR_k * relX;
					epsAdd.ref_r()[i] = epsMax - fR_eps * relX;
				}
			}

			epsAdd.ref_r() = std::pow(kAdd.ref(), 3. / 2.) / (lBound + rBound);

			phase->kTurb.ref_r() += kAdd.ref();
			phase->epsTurb.ref_r() += epsAdd.ref();

			std::cout << "Added linear peak for k and epsilon." << std::endl;
		}
		else if (profileType != "no")
			throw exception(
					"Wrong parameter of peak, must be <<linear>> or <<no>>.",
					errorsEnum::initializationError);
	}
	/**/

	/*Check for minimum values.*/
	std::replace_if(std::begin(phase->kTurb.ref_r()),
			std::end(phase->kTurb.ref_r()), [&phase](const scalar value) 
			{
				return value < phase->turbulenceSources->turbPar->mink();
			}, phase->turbulenceSources->turbPar->mink());

	std::replace_if(std::begin(phase->epsTurb.ref_r()),
			std::end(phase->epsTurb.ref_r()), [&phase](const scalar value) 
			{
				return value < phase->turbulenceSources->turbPar->mineps();
			}, phase->turbulenceSources->turbPar->mineps());

	/*Calculate derived fields.*/
	for (std::size_t k = 1; k < phase->concentration.v.size(); ++k)
		phase->concentration.v[0].ref_r() += phase->concentration.v[k].ref();

	for (std::size_t k = 1; k < phase->density.size(); ++k)
		phase->density[k].ref_r() = phase->concentration.v[k].ref()
				* phase->phaseThermodynamics->Mv()[k - 1];

	for (std::size_t k = 1; k < phase->density.size(); ++k)
		phase->density[0].ref_r() += phase->density[k].ref();

	phase->momentum.ref_r() =
			astProduct(phase->velocity, phase->density[0]).ref();

	phase->internalEnergy.ref_r() = phase->phaseThermodynamics->UvFromp(
			phase->concentration.p, phase->pressure.ref());

	phase->temperature.ref_r() = phase->phaseThermodynamics->TFromUv(
			phase->concentration.p, phase->internalEnergy.ref());

	phase->rhokTurb.ref_r() = phase->kTurb.ref() * phase->density[0].ref();

	phase->rhoepsTurb.ref_r() = phase->density[0].ref() * phase->epsTurb.ref();

	phase->rhoaTurb.ref_r() = astProduct(phase->aTurb, phase->density[0]).ref();

	phase->rhobTurb.ref_r() = phase->density[0].ref() * phase->bTurb.ref();

	phase->recalculateCoefficients(phase->kTurb, phase->epsTurb,
			*phase->turbulenceSources->turbPar);

	/*Set turbulent quantities to zero, if it is non-turbulent task.*/
	if (!phase->turbulenceSources->turbulence)
	{
		phase->kTurb.ref_r() = scalar { 0 };
		phase->epsTurb.ref_r() = 0;
		phase->rhokTurb.ref_r() = scalar { 0 };
		phase->rhoepsTurb.ref_r() = 0;
		phase->aTurb.ref_r() = vector(0);
		phase->bTurb.ref_r() = 0;
		phase->rhoaTurb.ref_r() = vector(0);
		phase->rhobTurb.ref_r() = 0;

		phase->tAssign(0.);
	}

	{
		const auto v2 = ampProduct(phase->velocity, phase->velocity);

		phase->totalEnergy.ref_r() = phase->internalEnergy.ref()
				+ phase->density[0].ref() * v2.ref() * 0.5
				+ phase->rhokTurb.ref();
	}

	phase->HelmholtzEnergy.ref_r() = phase->phaseThermodynamics->Fv(
			phase->concentration.p, phase->temperature.ref());

	phase->entropy.ref_r() = phase->phaseThermodynamics->Sv(
			phase->concentration.p, phase->temperature.ref());

	return std::make_tuple(std::move(phase), enthalpyFlowFlag,
			molMassDiffusionFlag);
}
