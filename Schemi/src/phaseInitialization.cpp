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

	std::vector<std::vector<std::string>> matrixOfSubstancesConditions(
			numberOfComponents,
			std::vector<std::string>(4 + 12 + numberOfZones));

	transportCoefficients<cubicCell> tCoeffsPhase { meshReference };

	bunchOfFields<cubicCell> cellFields { meshReference, numberOfComponents };
	std::unique_ptr<abstractMixtureThermodynamics> mixture;

	/*Set boundary condition for transport coefficients.*/
	auto bndConCom = commonConditions;

	tCoeffsPhase.tD = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);
	tCoeffsPhase.tKappa = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);
	tCoeffsPhase.tLambda = volumeField<scalar>(meshReference, 0, bndConCom[0],
			0, bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4],
			0, bndConCom[5], 0);
	tCoeffsPhase.k_D = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);
	tCoeffsPhase.eps_D = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);
	tCoeffsPhase.a_D = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);
	tCoeffsPhase.b_D = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);

	std::replace(bndConCom.begin(), bndConCom.end(),
			boundaryConditionType::calculated,
			boundaryConditionType::calculatedTurbulentViscosity);
	tCoeffsPhase.tNu = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);

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
	tCoeffsPhase.pNu = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);
	tCoeffsPhase.pD = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);
	tCoeffsPhase.pKappa = volumeField<scalar>(meshReference, 0, bndConCom[0], 0,
			bndConCom[1], 0, bndConCom[2], 0, bndConCom[3], 0, bndConCom[4], 0,
			bndConCom[5], 0);
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
				>> skipBuffer >> matrixOfSubstancesConditions[k][3] /*Pcrit*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][4] /*Boundary condition type: tail*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][5] /*Boundary condition value: tail*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][6] /*Boundary condition type: point*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][7] /*Boundary condition value: point*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][8] /*Boundary condition type: bottom*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][9] /*Boundary condition value: bottom*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][10] /*Boundary condition type: right*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][11] /*Boundary condition value: right*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][12] /*Boundary condition type: left*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][13] /*Boundary condition value: left*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][14] /*Boundary condition type: top*/
				>> skipBuffer >> matrixOfSubstancesConditions[k][15]; /*Boundary condition value: top*/

		substanceConditionsFile >> skipBuffer;
		for (std::size_t j = 0; j < numberOfZones; ++j)
			substanceConditionsFile >> matrixOfSubstancesConditions[k][16 + j]; /*Values in zones*/
		substanceConditionsFile.close();
	}

	/*Convert strings to values.*/
	std::valarray<std::valarray<scalar>> thermodynamicalProperties {
			std::valarray<scalar>(numberOfComponents), 4 };

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
		std::vector<boundaryConditionType> bndCon(6);
		std::vector<std::vector<scalar>> bndVCon(6);

		for (std::size_t j = 0; j < bndCon.size(); ++j)
			boundaryConditionFromString(
					matrixOfSubstancesConditions[k - 1][4 + 2 * j],
					matrixOfSubstancesConditions[k - 1][5 + 2 * j], bndCon[j],
					bndVCon[j]);

		parallelism.correctBoundaryConditions(bndCon);

		cellFields.concentration.v[k] = volumeField<scalar>(meshReference, 0,
				bndCon[0], bndVCon[0], bndCon[1], bndVCon[1], bndCon[2],
				bndVCon[2], bndCon[3], bndVCon[3], bndCon[4], bndVCon[4],
				bndCon[5], bndVCon[5], parallelism.mpi_size);

		auto bndVConRho_k = bndVCon;

		for (auto & b_s : bndVConRho_k)
			for (auto & b_i : b_s)
				b_i *= thermodynamicalProperties[0][k - 1];

		cellFields.density[k] = volumeField<scalar>(meshReference, 0, bndCon[0],
				bndVConRho_k[0], bndCon[1], bndVConRho_k[1], bndCon[2],
				bndVConRho_k[2], bndCon[3], bndVConRho_k[3], bndCon[4],
				bndVConRho_k[4], bndCon[5], bndVConRho_k[5],
				parallelism.mpi_size);
	}
	{
		cellFields.concentration.v[0] = volumeField<scalar>(meshReference, 0);

		auto bndConRho = commonConditions;

		std::replace(bndConRho.begin(), bndConRho.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedDensity);

		cellFields.density[0] = volumeField<scalar>(meshReference, 0,
				bndConRho[0], 0, bndConRho[1], 0, bndConRho[2], 0, bndConRho[3],
				0, bndConRho[4], 0, bndConRho[5], 0);
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

	std::vector<std::string> velocityConditions(24 + 3 * numberOfZones);
	{
		std::ifstream velocityConditionsFile { "./set/velocity.txt",
				std::ios::in };

		if (velocityConditionsFile.is_open())
			std::cout << "./set/velocity.txt is opened." << std::endl;
		else
			throw exception("./set/velocity.txt not found.",
					errorsEnum::initializationError);

		velocityConditionsFile >> skipBuffer >> velocityConditions[0]
				>> skipBuffer >> velocityConditions[1] >> velocityConditions[2]
				>> velocityConditions[3] >> skipBuffer >> velocityConditions[4]
				>> skipBuffer >> velocityConditions[5] >> velocityConditions[6]
				>> velocityConditions[7] >> skipBuffer >> velocityConditions[8]
				>> skipBuffer >> velocityConditions[9] >> velocityConditions[10]
				>> velocityConditions[11] >> skipBuffer
				>> velocityConditions[12] >> skipBuffer
				>> velocityConditions[13] >> velocityConditions[14]
				>> velocityConditions[15] >> skipBuffer
				>> velocityConditions[16] >> skipBuffer
				>> velocityConditions[17] >> velocityConditions[18]
				>> velocityConditions[19] >> skipBuffer
				>> velocityConditions[20] >> skipBuffer
				>> velocityConditions[21] >> velocityConditions[22]
				>> velocityConditions[23];

		velocityConditionsFile >> skipBuffer;
		for (std::size_t i = 0; i < numberOfZones; ++i)
			velocityConditionsFile >> velocityConditions[24 + 3 * i]
					>> velocityConditions[25 + 3 * i]
					>> velocityConditions[26 + 3 * i];

		velocityConditionsFile.close();

		std::vector<boundaryConditionType> bndCon(6);
		std::vector<std::vector<vector>> bndVCon(6);

		for (std::size_t j = 0; j < bndCon.size(); ++j)
			boundaryConditionFromString(velocityConditions[0 + 4 * j],
					velocityConditions[1 + 4 * j],
					velocityConditions[2 + 4 * j],
					velocityConditions[3 + 4 * j], bndCon[j], bndVCon[j]);

		parallelism.correctBoundaryConditions(bndCon);

		phase->velocity = volumeField<vector>(meshReference, vector(0),
				bndCon[0], bndVCon[0], bndCon[1], bndVCon[1], bndCon[2],
				bndVCon[2], bndCon[3], bndVCon[3], bndCon[4], bndVCon[4],
				bndCon[5], bndVCon[5], parallelism.mpi_size);
	}

	std::vector<std::string> pressureConditions(12 + numberOfZones);
	{
		std::ifstream pressureConditionsFile { "./set/pressure.txt",
				std::ios::in };

		if (pressureConditionsFile.is_open())
			std::cout << "./set/pressure.txt is opened." << std::endl;
		else
			throw exception("./set/pressure.txt not found.",
					errorsEnum::initializationError);

		pressureConditionsFile >> skipBuffer >> pressureConditions[0]
				>> skipBuffer >> pressureConditions[1] >> skipBuffer
				>> pressureConditions[2] >> skipBuffer >> pressureConditions[3]
				>> skipBuffer >> pressureConditions[4] >> skipBuffer
				>> pressureConditions[5] >> skipBuffer >> pressureConditions[6]
				>> skipBuffer >> pressureConditions[7] >> skipBuffer
				>> pressureConditions[8] >> skipBuffer >> pressureConditions[9]
				>> skipBuffer >> pressureConditions[10] >> skipBuffer
				>> pressureConditions[11];

		pressureConditionsFile >> skipBuffer;
		for (std::size_t i = 0; i < numberOfZones; ++i)
			pressureConditionsFile >> pressureConditions[12 + i];

		pressureConditionsFile.close();

		std::vector<boundaryConditionType> bndCon(6);
		std::vector<std::vector<scalar>> bndVCon(6);

		for (std::size_t j = 0; j < bndCon.size(); ++j)
			boundaryConditionFromString(pressureConditions[0 + 2 * j],
					pressureConditions[1 + 2 * j], bndCon[j], bndVCon[j]);

		parallelism.correctBoundaryConditions(bndCon);

		phase->pressure = volumeField<scalar>(meshReference, 0, bndCon[0],
				bndVCon[0], bndCon[1], bndVCon[1], bndCon[2], bndVCon[2],
				bndCon[3], bndVCon[3], bndCon[4], bndVCon[4], bndCon[5],
				bndVCon[5], parallelism.mpi_size);
	}

	std::vector<std::string> kTurbConditions(12 + numberOfZones);
	{
		std::ifstream kTurbConditionsFile { "./set/k.txt", std::ios::in };

		if (kTurbConditionsFile.is_open())
			std::cout << "./set/k.txt is opened." << std::endl;
		else
			throw exception("./set/k.txt is opened.",
					errorsEnum::initializationError);

		kTurbConditionsFile >> skipBuffer >> kTurbConditions[0] >> skipBuffer
				>> kTurbConditions[1] >> skipBuffer >> kTurbConditions[2]
				>> skipBuffer >> kTurbConditions[3] >> skipBuffer
				>> kTurbConditions[4] >> skipBuffer >> kTurbConditions[5]
				>> skipBuffer >> kTurbConditions[6] >> skipBuffer
				>> kTurbConditions[7] >> skipBuffer >> kTurbConditions[8]
				>> skipBuffer >> kTurbConditions[9] >> skipBuffer
				>> kTurbConditions[10] >> skipBuffer >> kTurbConditions[11];

		kTurbConditionsFile >> skipBuffer;
		for (std::size_t i = 0; i < numberOfZones; ++i)
			kTurbConditionsFile >> kTurbConditions[12 + i];

		kTurbConditionsFile.close();

		std::vector<boundaryConditionType> bndCon(6);
		std::vector<std::vector<scalar>> bndVCon(6);

		for (std::size_t j = 0; j < bndCon.size(); ++j)
			boundaryConditionFromString(kTurbConditions[0 + 2 * j],
					kTurbConditions[1 + 2 * j], bndCon[j], bndVCon[j]);

		parallelism.correctBoundaryConditions(bndCon);

		phase->kTurb = volumeField<scalar>(meshReference, scalar { 0 },
				bndCon[0], bndVCon[0], bndCon[1], bndVCon[1], bndCon[2],
				bndVCon[2], bndCon[3], bndVCon[3], bndCon[4], bndVCon[4],
				bndCon[5], bndVCon[5], parallelism.mpi_size);
	}

	std::vector<std::string> epsTurbConditions(12 + numberOfZones);
	{
		std::ifstream epsTurbConditionsFile { "./set/epsilon.txt", std::ios::in };

		if (epsTurbConditionsFile.is_open())
			std::cout << "./set/epsilon.txt is opened." << std::endl;
		else
			throw exception("./set/epsilon.txt not found.",
					errorsEnum::initializationError);

		epsTurbConditionsFile >> skipBuffer >> epsTurbConditions[0]
				>> skipBuffer >> epsTurbConditions[1] >> skipBuffer
				>> epsTurbConditions[2] >> skipBuffer >> epsTurbConditions[3]
				>> skipBuffer >> epsTurbConditions[4] >> skipBuffer
				>> epsTurbConditions[5] >> skipBuffer >> epsTurbConditions[6]
				>> skipBuffer >> epsTurbConditions[7] >> skipBuffer
				>> epsTurbConditions[8] >> skipBuffer >> epsTurbConditions[9]
				>> skipBuffer >> epsTurbConditions[10] >> skipBuffer
				>> epsTurbConditions[11];

		epsTurbConditionsFile >> skipBuffer;
		for (std::size_t i = 0; i < numberOfZones; ++i)
			epsTurbConditionsFile >> epsTurbConditions[12 + i];

		epsTurbConditionsFile.close();

		std::vector<boundaryConditionType> bndCon(6);
		std::vector<std::vector<scalar>> bndVCon(6);

		for (std::size_t j = 0; j < bndCon.size(); ++j)
			boundaryConditionFromString(epsTurbConditions[0 + 2 * j],
					epsTurbConditions[1 + 2 * j], bndCon[j], bndVCon[j]);

		parallelism.correctBoundaryConditions(bndCon);

		phase->epsTurb = volumeField<scalar>(meshReference, 0, bndCon[0],
				bndVCon[0], bndCon[1], bndVCon[1], bndCon[2], bndVCon[2],
				bndCon[3], bndVCon[3], bndCon[4], bndVCon[4], bndCon[5],
				bndVCon[5], parallelism.mpi_size);
	}

	std::vector<std::string> aTurbConditions(24 + 3 * numberOfZones);
	{
		std::ifstream aTurbConditionsFile { "./set/a.txt", std::ios::in };

		if (aTurbConditionsFile.is_open())
			std::cout << "./set/a.txt is opened." << std::endl;
		else
			throw exception("./set/a.txt not found.",
					errorsEnum::initializationError);

		aTurbConditionsFile >> skipBuffer >> aTurbConditions[0] >> skipBuffer
				>> aTurbConditions[1] >> aTurbConditions[2]
				>> aTurbConditions[3] >> skipBuffer >> aTurbConditions[4]
				>> skipBuffer >> aTurbConditions[5] >> aTurbConditions[6]
				>> aTurbConditions[7] >> skipBuffer >> aTurbConditions[8]
				>> skipBuffer >> aTurbConditions[9] >> aTurbConditions[10]
				>> aTurbConditions[11] >> skipBuffer >> aTurbConditions[12]
				>> skipBuffer >> aTurbConditions[13] >> aTurbConditions[14]
				>> aTurbConditions[15] >> skipBuffer >> aTurbConditions[16]
				>> skipBuffer >> aTurbConditions[17] >> aTurbConditions[18]
				>> aTurbConditions[19] >> skipBuffer >> aTurbConditions[20]
				>> skipBuffer >> aTurbConditions[21] >> aTurbConditions[22]
				>> aTurbConditions[23];

		aTurbConditionsFile >> skipBuffer;
		for (std::size_t i = 0; i < numberOfZones; ++i)
			aTurbConditionsFile >> aTurbConditions[24 + 3 * i]
					>> aTurbConditions[25 + 3 * i]
					>> aTurbConditions[26 + 3 * i];

		aTurbConditionsFile.close();

		std::vector<boundaryConditionType> bndCon(6);
		std::vector<std::vector<vector>> bndVCon(6);

		for (std::size_t j = 0; j < bndCon.size(); ++j)
			boundaryConditionFromString(aTurbConditions[0 + 4 * j],
					aTurbConditions[1 + 2 * j], aTurbConditions[2 + 2 * j],
					aTurbConditions[3 + 2 * j], bndCon[j], bndVCon[j]);

		parallelism.correctBoundaryConditions(bndCon);

		phase->aTurb = volumeField<vector>(meshReference, vector(0, 0, 0),
				bndCon[0], bndVCon[0], bndCon[1], bndVCon[1], bndCon[2],
				bndVCon[2], bndCon[3], bndVCon[3], bndCon[4], bndVCon[4],
				bndCon[5], bndVCon[5], parallelism.mpi_size);
	}

	std::vector<std::string> bTurbConditions(12 + numberOfZones);
	{
		std::ifstream bTurbConditionsFile { "./set/b.txt", std::ios::in };

		if (bTurbConditionsFile.is_open())
			std::cout << "./set/b.txt is opened." << std::endl;
		else
			throw exception("./set/b.txt not found.",
					errorsEnum::initializationError);

		bTurbConditionsFile >> skipBuffer >> bTurbConditions[0] >> skipBuffer
				>> bTurbConditions[1] >> skipBuffer >> bTurbConditions[2]
				>> skipBuffer >> bTurbConditions[3] >> skipBuffer
				>> bTurbConditions[4] >> skipBuffer >> bTurbConditions[5]
				>> skipBuffer >> bTurbConditions[6] >> skipBuffer
				>> bTurbConditions[7] >> skipBuffer >> bTurbConditions[8]
				>> skipBuffer >> bTurbConditions[9] >> skipBuffer
				>> bTurbConditions[10] >> skipBuffer >> bTurbConditions[11];

		bTurbConditionsFile >> skipBuffer;
		for (std::size_t i = 0; i < numberOfZones; ++i)
			bTurbConditionsFile >> bTurbConditions[12 + i];

		bTurbConditionsFile.close();

		std::vector<boundaryConditionType> bndCon(6);
		std::vector<std::vector<scalar>> bndVCon(6);

		for (std::size_t j = 0; j < bndCon.size(); ++j)
			boundaryConditionFromString(bTurbConditions[0 + 2 * j],
					bTurbConditions[1 + 2 * j], bndCon[j], bndVCon[j]);

		parallelism.correctBoundaryConditions(bndCon);

		phase->bTurb = volumeField<scalar>(meshReference,
				phase->turbulenceSources->turbPar->minb_value, bndCon[0],
				bndVCon[0], bndCon[1], bndVCon[1], bndCon[2], bndVCon[2],
				bndCon[3], bndVCon[3], bndCon[4], bndVCon[4], bndCon[5],
				bndVCon[5], parallelism.mpi_size);
	}

	/*Creating other fields.*/
	phase->momentum = volumeField<vector>(meshReference, vector(0),
			commonConditions[0], vector(0), commonConditions[1], vector(0),
			commonConditions[2], vector(0), commonConditions[3], vector(0),
			commonConditions[4], vector(0), commonConditions[5], vector(0));

	phase->internalEnergy = volumeField<scalar>(meshReference, 0,
			commonConditions[0], 0, commonConditions[1], 0, commonConditions[2],
			0, commonConditions[3], 0, commonConditions[4], 0,
			commonConditions[5], 0);

	phase->totalEnergy = volumeField<scalar>(meshReference, 0,
			commonConditions[0], 0, commonConditions[1], 0, commonConditions[2],
			0, commonConditions[3], 0, commonConditions[4], 0,
			commonConditions[5], 0);

	{
		auto bndCon = commonConditions;
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedTemperature);
		phase->temperature = volumeField<scalar>(meshReference, 0, bndCon[0], 0,
				bndCon[1], 0, bndCon[2], 0, bndCon[3], 0, bndCon[4], 0,
				bndCon[5], 0);
	}

	phase->rhokTurb = volumeField<scalar>(meshReference, scalar { 0 },
			commonConditions[0], scalar { 0 }, commonConditions[1],
			scalar { 0 }, commonConditions[2], scalar { 0 },
			commonConditions[3], scalar { 0 }, commonConditions[4],
			scalar { 0 }, commonConditions[5], scalar { 0 });

	phase->rhoepsTurb = volumeField<scalar>(meshReference, 0,
			commonConditions[0], 0, commonConditions[1], 0, commonConditions[2],
			0, commonConditions[3], 0, commonConditions[4], 0,
			commonConditions[5], 0);

	phase->rhoaTurb = volumeField<vector>(meshReference, vector(0),
			commonConditions[0], vector(0), commonConditions[1], vector(0),
			commonConditions[2], vector(0), commonConditions[3], vector(0),
			commonConditions[4], vector(0), commonConditions[5], vector(0));

	phase->rhobTurb = volumeField<scalar>(meshReference, 0, commonConditions[0],
			0, commonConditions[1], 0, commonConditions[2], 0,
			commonConditions[3], 0, commonConditions[4], 0, commonConditions[5],
			0);

	phase->HelmholtzEnergy = volumeField<scalar>(meshReference, 0,
			commonConditions[0], 0, commonConditions[1], 0, commonConditions[2],
			0, commonConditions[3], 0, commonConditions[4], 0,
			commonConditions[5], 0);

	phase->entropy = volumeField<scalar>(meshReference, 0, commonConditions[0],
			0, commonConditions[1], 0, commonConditions[2], 0,
			commonConditions[3], 0, commonConditions[4], 0, commonConditions[5],
			0);

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
				if ((radiusOfCell.v()[0] >= vectorsOfZones[2 * j].v()[0])
						&& (radiusOfCell.v()[0]
								< vectorsOfZones[2 * j + 1].v()[0])
						&& (radiusOfCell.v()[1] >= vectorsOfZones[2 * j].v()[1])
						&& (radiusOfCell.v()[1]
								< vectorsOfZones[2 * j + 1].v()[1])
						&& (radiusOfCell.v()[2] >= vectorsOfZones[2 * j].v()[2])
						&& (radiusOfCell.v()[2]
								<= vectorsOfZones[2 * j + 1].v()[2]))
				{
					zoneFounded = true;
					phase->concentration.v[k].ref_r()[i] = std::stod(
							matrixOfSubstancesConditions[k - 1][16 + j]);
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
			if ((radiusOfCell.v()[0] >= vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] >= vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] >= vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->velocity.ref_r()[i] = vector(
						std::stod(velocityConditions[24 + 3 * j]),
						std::stod(velocityConditions[25 + 3 * j]),
						std::stod(velocityConditions[26 + 3 * j]));
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
			if ((radiusOfCell.v()[0] >= vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] >= vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] >= vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->pressure.ref_r()[i] = std::stod(
						pressureConditions[12 + j]);
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
			if ((radiusOfCell.v()[0] >= vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] >= vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] >= vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->kTurb.ref_r()[i] = scalar { std::stod(
						kTurbConditions[12 + j]) };
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
			if ((radiusOfCell.v()[0] >= vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] >= vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] >= vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->epsTurb.ref_r()[i] = std::stod(
						epsTurbConditions[12 + j]);
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
			if ((radiusOfCell.v()[0] >= vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] >= vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] >= vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->aTurb.ref_r()[i] = vector(
						std::stod(aTurbConditions[24 + 3 * j]),
						std::stod(aTurbConditions[25 + 3 * j]),
						std::stod(aTurbConditions[26 + 3 * j]));
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
			if ((radiusOfCell.v()[0] >= vectorsOfZones[2 * j].v()[0])
					&& (radiusOfCell.v()[0] < vectorsOfZones[2 * j + 1].v()[0])
					&& (radiusOfCell.v()[1] >= vectorsOfZones[2 * j].v()[1])
					&& (radiusOfCell.v()[1] < vectorsOfZones[2 * j + 1].v()[1])
					&& (radiusOfCell.v()[2] >= vectorsOfZones[2 * j].v()[2])
					&& (radiusOfCell.v()[2] < vectorsOfZones[2 * j + 1].v()[2]))
			{
				zoneFounded = true;
				phase->bTurb.ref_r()[i] = std::stod(bTurbConditions[12 + j]);
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
