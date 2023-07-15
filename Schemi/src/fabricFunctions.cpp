/*
 * fabricFunctions.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "fabricFunctions.hpp"

#include <iostream>

#include "typeOfTVDLimiterEnum.hpp"
#include "flowSolverEnum.hpp"

#include "arithmeticAGen.hpp"
#include "BHRGen.hpp"
#include "BHRKLGen.hpp"
#include "chemicalKineticsChlorumDissociation.hpp"
#include "chemicalKineticsChlorumHydrogeniumDissociation.hpp"
#include "chemicalKineticsH2Cl2Combustion.hpp"
#include "chemicalKineticsNoReaction.hpp"
#include "chemicalKineticsNO2Disproportionation.hpp"
#include "conjugateGradientSolver.hpp"
#include "decayGen.hpp"
#include "GaussSeidelSolver.hpp"
#include "HLLC2pSolver.hpp"
#include "HLLCFSolver.hpp"
#include "HLLCLMSolver.hpp"
#include "HLLCSolver.hpp"
#include "HLLSolver.hpp"
#include "HQUICKLimiter.hpp"
#include "JacobiConjugateGradientSolver.hpp"
#include "kEpsAGen.hpp"
#include "KTSolver.hpp"
#include "linearLimiter.hpp"
#include "minmod2Limiter.hpp"
#include "minmodLimiter.hpp"
#include "noMatrixSolver.hpp"
#include "RichtmyerSolver.hpp"
#include "shearGen.hpp"
#include "superbeeLimiter.hpp"
#include "ThomasSolver.hpp"
#include "vanAlbada2Limiter.hpp"
#include "vanAlbadaLimiter.hpp"
#include "vanLeer2Limiter.hpp"
#include "vanLeerLimiter.hpp"
#include "zeroGen.hpp"
#include "zeroLimiter.hpp"

std::unique_ptr<schemi::abstractLimiter> schemi::createLimiter(
		const std::string & name)
{
	typeOfTVDLimiterEnum limiterFlag;
	if (name == "minmod")
		limiterFlag = typeOfTVDLimiterEnum::minmod;
	else if (name == "vanLeer")
		limiterFlag = typeOfTVDLimiterEnum::vanLeer;
	else if (name == "linear")
		limiterFlag = typeOfTVDLimiterEnum::linear;
	else if (name == "zero")
		limiterFlag = typeOfTVDLimiterEnum::zero;
	else if (name == "vanAlbada")
		limiterFlag = typeOfTVDLimiterEnum::vanAlbada;
	else if (name == "HQUICK")
		limiterFlag = typeOfTVDLimiterEnum::HQUICK;
	else if (name == "vanLeer2")
		limiterFlag = typeOfTVDLimiterEnum::vanLeer2;
	else if (name == "superbee")
		limiterFlag = typeOfTVDLimiterEnum::superbee;
	else if (name == "vanAlbada2")
		limiterFlag = typeOfTVDLimiterEnum::vanAlbada2;
	else if (name == "minmod2")
		limiterFlag = typeOfTVDLimiterEnum::minmod2;
	else
		throw exception("Unknown TVD limiter flag.",
				errorsEnum::initializationError);

	switch (limiterFlag)
	{
	case typeOfTVDLimiterEnum::minmod:
		return std::make_unique<minmodLimiter>();
		break;
	case typeOfTVDLimiterEnum::vanLeer:
		return std::make_unique<vanLeerLimiter>();
		break;
	case typeOfTVDLimiterEnum::linear:
		return std::make_unique<linearLimiter>();
		break;
	case typeOfTVDLimiterEnum::vanAlbada:
		return std::make_unique<vanAlbadaLimiter>();
		break;
	case typeOfTVDLimiterEnum::HQUICK:
		return std::make_unique<HQUICKLimiter>();
		break;
	case typeOfTVDLimiterEnum::vanLeer2:
		return std::make_unique<vanLeer2Limiter>();
		break;
	case typeOfTVDLimiterEnum::superbee:
		return std::make_unique<superbeeLimiter>();
		break;
	case typeOfTVDLimiterEnum::vanAlbada2:
		return std::make_unique<vanAlbada2Limiter>();
		break;
	case typeOfTVDLimiterEnum::minmod2:
		return std::make_unique<minmod2Limiter>();
		break;
	default:
		return std::make_unique<zeroLimiter>();
		break;
	}
}

std::pair<std::unique_ptr<schemi::abstractMatrixSolver>,
		std::unique_ptr<schemi::abstractMatrixSolver>> schemi::createMatrixSolver(
		const std::string & name, const std::string & dim,
		const std::size_t iter)
{
	matrixSolverEnum matrixSolverFlag;
	if (name == "GaussSeidel")
		matrixSolverFlag = matrixSolverEnum::GaussSeidel;
	else if ((name == "Thomas") && (dim == "1D"))
		matrixSolverFlag = matrixSolverEnum::Thomas;
	else if (name == "noSolver")
		matrixSolverFlag = matrixSolverEnum::noSovler;
	else if (name == "explicit")
		matrixSolverFlag = matrixSolverEnum::explicitSolver;
	else if (name == "conjugateGradient")
		matrixSolverFlag = matrixSolverEnum::conjugateGradient;
	else if (name == "JacobiConjugateGradient")
		matrixSolverFlag = matrixSolverEnum::JacobiConjugateGradient;
	else
		throw exception("Unknown or inappropriate matrix solver flag.",
				errorsEnum::initializationError);

	switch (matrixSolverFlag)
	{
	case matrixSolverEnum::conjugateGradient:
		return std::make_pair(
				std::make_unique<conjugateGradientSovler>(iter,
						matrixSolverFlag),
				std::make_unique<conjugateGradientSovler>(iter,
						matrixSolverFlag));
		break;
	case matrixSolverEnum::JacobiConjugateGradient:
		return std::make_pair(
				std::make_unique<JacobiConjugateGradientSolver>(iter,
						matrixSolverFlag),
				std::make_unique<JacobiConjugateGradientSolver>(iter,
						matrixSolverFlag));
		break;
	case matrixSolverEnum::Thomas:
		return std::make_pair(std::make_unique<ThomasSolver>(matrixSolverFlag),
				std::make_unique<GaussSeidelSolver>(iter,
						matrixSolverEnum::GaussSeidel));
		break;
	case matrixSolverEnum::explicitSolver:
	case matrixSolverEnum::noSovler:
		return std::make_pair(
				std::make_unique<noMatrixSolver>(matrixSolverFlag),
				std::make_unique<noMatrixSolver>(matrixSolverFlag));
		break;
	default:
		return std::make_pair(
				std::make_unique<GaussSeidelSolver>(iter,
						matrixSolverEnum::GaussSeidel),
				std::make_unique<GaussSeidelSolver>(iter,
						matrixSolverEnum::GaussSeidel));
		break;
	}
}

std::unique_ptr<schemi::abstractFlowSolver> schemi::createFlowSolver(
		const std::string & name)
{
	flowSolverEnum flowSolwerFlag;
	if (name == "HLLCF")
		flowSolwerFlag = flowSolverEnum::HLLCF;
	else if (name == "HLLC")
		flowSolwerFlag = flowSolverEnum::HLLC;
	else if (name == "HLL")
		flowSolwerFlag = flowSolverEnum::HLL;
	else if (name == "HLLCLM")
		flowSolwerFlag = flowSolverEnum::HLLCLM;
	else if (name == "KT")
		flowSolwerFlag = flowSolverEnum::KT;
	else if (name == "HLLC2p")
		flowSolwerFlag = flowSolverEnum::HLLC2p;
	else if (name == "Richtmyer")
		flowSolwerFlag = flowSolverEnum::Richtmyer;
	else
		throw exception("Unknown flow solver flag.",
				errorsEnum::initializationError);

	switch (flowSolwerFlag)
	{
	case flowSolverEnum::HLL:
		return std::make_unique<HLLSolver>();
		break;
	case flowSolverEnum::HLLCF:
		return std::make_unique<HLLCFSolver>();
		break;
	case flowSolverEnum::HLLC:
		return std::make_unique<HLLCSolver>();
		break;
	case flowSolverEnum::HLLCLM:
		return std::make_unique<HLLCLMSolver>();
		break;
	case flowSolverEnum::HLLC2p:
		return std::make_unique<HLLC2pSolver>();
		break;
	case flowSolverEnum::Richtmyer:
		return std::make_unique<RichtmyerSolver>();
		break;
	default:
		return std::make_unique<KTSolver>();
		break;
	}
}

std::unique_ptr<schemi::abstractTurbulenceGen> schemi::createTurbulenceModel(
		const std::string & turbulenceONString,
		const std::string & sourceTypeString)
{
	bool turbulenceFlag;
	if (turbulenceONString == "on")
		turbulenceFlag = true;
	else if (turbulenceONString == "off")
		turbulenceFlag = false;
	else
		throw exception("Unknown turbulence flag.",
				errorsEnum::initializationError);

	turbulenceModelEnum turbulenceModelFlag;
	if (sourceTypeString == "BHR")
		turbulenceModelFlag = turbulenceModelEnum::BHRSource;
	else if (sourceTypeString == "zero")
		turbulenceModelFlag = turbulenceModelEnum::zeroSource;
	else if (sourceTypeString == "decay")
		turbulenceModelFlag = turbulenceModelEnum::decaySource;
	else if (sourceTypeString == "shear")
		turbulenceModelFlag = turbulenceModelEnum::shearSource;
	else if (sourceTypeString == "arithmetic1")
		turbulenceModelFlag = turbulenceModelEnum::arithmeticA1Source;
	else if (sourceTypeString == "arithmetic2")
		turbulenceModelFlag = turbulenceModelEnum::arithmeticA2Source;
	else if (sourceTypeString == "kEpsA")
		turbulenceModelFlag = turbulenceModelEnum::kEpsASource;
	else if (sourceTypeString == "arithmetic3")
		turbulenceModelFlag = turbulenceModelEnum::arithmeticA3Source;
	else if (sourceTypeString == "BHRKL")
		turbulenceModelFlag = turbulenceModelEnum::BHRKLSource;
	else
		throw exception("Unknown source generation flag.",
				errorsEnum::initializationError);

	std::unique_ptr<abstractTurbulenceGen> turbulenceSources;
	switch (turbulenceModelFlag)
	{
	case turbulenceModelEnum::BHRSource:
		turbulenceSources = std::make_unique<BHRGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModelEnum::decaySource:
		turbulenceSources = std::make_unique<decayGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModelEnum::shearSource:
		turbulenceSources = std::make_unique<shearGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModelEnum::arithmeticA1Source:
	case turbulenceModelEnum::arithmeticA2Source:
	case turbulenceModelEnum::arithmeticA3Source:
		turbulenceSources = std::make_unique<arithmeticAGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModelEnum::kEpsASource:
		turbulenceSources = std::make_unique<kEpsAGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModelEnum::BHRKLSource:
		turbulenceSources = std::make_unique<BHRKLGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	default:
		turbulenceSources = std::make_unique<zeroGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	}

	/*Read turbulent parameters.*/
	{
		std::ifstream turbulentParametersFile { "./set/turbulentParameters.txt" };
		if (turbulentParametersFile.is_open())
			std::cout << "./set/turbulentParameters.txt is opened."
					<< std::endl;
		else
			throw exception("./set/turbulentParameters.txt not found.",
					errorsEnum::initializationError);

		turbulenceSources->turbPar->readTurbulentParameters(
				turbulentParametersFile);

		turbulentParametersFile.close();
	}

	return turbulenceSources;
}

std::unique_ptr<schemi::abstractChemicalKinetics> schemi::createChemicalKinetics(
		const homogeneousPhase<cubicCell> & phaseIn,
		const chemicalReactionsEnum chemReactFlag)
{
	switch (chemReactFlag)
	{
	case chemicalReactionsEnum::Cl2Dissociation:
		return std::make_unique<chemicalKineticsChlorumDissociation>(phaseIn);
		break;
	case chemicalReactionsEnum::Cl2H2Dissociation:
		return std::make_unique<chemicalKineticsChlorumHydrogeniumDissociation>(
				phaseIn);
		break;
	case chemicalReactionsEnum::H2Cl2Combustion:
		return std::make_unique<chemicalKineticsH2Cl2Combustion>(phaseIn);
		break;
	case chemicalReactionsEnum::NO2Disproportionation:
		return std::make_unique<chemicalKineticsNO2Disproportionation>(phaseIn);
		break;
	case chemicalReactionsEnum::noReaction:
	default:
		return std::make_unique<chemicalKineticsNoReaction>();
		break;
	}
}
