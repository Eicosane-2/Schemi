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
#include "hardSpheresTransportModel.hpp"
#include "HLLC2pSolver.hpp"
#include "HLLCFSolver.hpp"
#include "HLLCLMSolver.hpp"
#include "HLLCSolver.hpp"
#include "HLLSolver.hpp"
#include "HQUICKLimiter.hpp"
#include "JacobiConjugateGradientSolver.hpp"
#include "JacobiSolver.hpp"
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
	typeOfTVDLimiter limiterFlag;
	if (name == "minmod")
		limiterFlag = typeOfTVDLimiter::minmod;
	else if (name == "vanLeer")
		limiterFlag = typeOfTVDLimiter::vanLeer;
	else if (name == "linear")
		limiterFlag = typeOfTVDLimiter::linear;
	else if (name == "zero")
		limiterFlag = typeOfTVDLimiter::zero;
	else if (name == "vanAlbada")
		limiterFlag = typeOfTVDLimiter::vanAlbada;
	else if (name == "HQUICK")
		limiterFlag = typeOfTVDLimiter::HQUICK;
	else if (name == "vanLeer2")
		limiterFlag = typeOfTVDLimiter::vanLeer2;
	else if (name == "superbee")
		limiterFlag = typeOfTVDLimiter::superbee;
	else if (name == "vanAlbada2")
		limiterFlag = typeOfTVDLimiter::vanAlbada2;
	else if (name == "minmod2")
		limiterFlag = typeOfTVDLimiter::minmod2;
	else
		throw exception("Unknown TVD limiter flag.",
				errors::initializationError);

	switch (limiterFlag)
	{
	case typeOfTVDLimiter::minmod:
		return std::make_unique<minmodLimiter>();
		break;
	case typeOfTVDLimiter::vanLeer:
		return std::make_unique<vanLeerLimiter>();
		break;
	case typeOfTVDLimiter::linear:
		return std::make_unique<linearLimiter>();
		break;
	case typeOfTVDLimiter::vanAlbada:
		return std::make_unique<vanAlbadaLimiter>();
		break;
	case typeOfTVDLimiter::HQUICK:
		return std::make_unique<HQUICKLimiter>();
		break;
	case typeOfTVDLimiter::vanLeer2:
		return std::make_unique<vanLeer2Limiter>();
		break;
	case typeOfTVDLimiter::superbee:
		return std::make_unique<superbeeLimiter>();
		break;
	case typeOfTVDLimiter::vanAlbada2:
		return std::make_unique<vanAlbada2Limiter>();
		break;
	case typeOfTVDLimiter::minmod2:
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
	matrixSolver matrixSolverFlag;
	if (name == "GaussSeidel")
		matrixSolverFlag = matrixSolver::GaussSeidel;
	else if ((name == "Thomas") && (dim == "1D"))
		matrixSolverFlag = matrixSolver::Thomas;
	else if (name == "noSolver")
		matrixSolverFlag = matrixSolver::noSovler;
	else if (name == "explicit")
		matrixSolverFlag = matrixSolver::explicitSolver;
	else if (name == "conjugateGradient")
		matrixSolverFlag = matrixSolver::conjugateGradient;
	else if (name == "JacobiConjugateGradient")
		matrixSolverFlag = matrixSolver::JacobiConjugateGradient;
	else if (name == "Jacobi")
		matrixSolverFlag = matrixSolver::Jacobi;
	else
		throw exception("Unknown or inappropriate matrix solver flag.",
				errors::initializationError);

	switch (matrixSolverFlag)
	{
	case matrixSolver::conjugateGradient:
		return std::make_pair(
				std::make_unique<conjugateGradientSovler>(iter,
						matrixSolverFlag),
				std::make_unique<conjugateGradientSovler>(iter,
						matrixSolverFlag));
		break;
	case matrixSolver::JacobiConjugateGradient:
		return std::make_pair(
				std::make_unique<JacobiConjugateGradientSolver>(iter,
						matrixSolverFlag),
				std::make_unique<JacobiConjugateGradientSolver>(iter,
						matrixSolverFlag));
		break;
	case matrixSolver::Thomas:
		return std::make_pair(std::make_unique<ThomasSolver>(matrixSolverFlag),
				std::make_unique<GaussSeidelSolver>(iter,
						matrixSolver::GaussSeidel));
		break;
	case matrixSolver::explicitSolver:
	case matrixSolver::noSovler:
		return std::make_pair(
				std::make_unique<noMatrixSolver>(matrixSolverFlag),
				std::make_unique<noMatrixSolver>(matrixSolverFlag));
		break;
	case matrixSolver::Jacobi:
		return std::make_pair(
				std::make_unique<JacobiSolver>(iter, matrixSolver::Jacobi),
				std::make_unique<JacobiSolver>(iter, matrixSolver::Jacobi));
		break;
	default:
		return std::make_pair(
				std::make_unique<GaussSeidelSolver>(iter,
						matrixSolver::GaussSeidel),
				std::make_unique<GaussSeidelSolver>(iter,
						matrixSolver::GaussSeidel));
		break;
	}
}

std::unique_ptr<schemi::abstractFlowSolver> schemi::createFlowSolver(
		const std::string & name)
{
	flowSolver flowSolwerFlag;
	if (name == "HLLCF")
		flowSolwerFlag = flowSolver::HLLCF;
	else if (name == "HLLC")
		flowSolwerFlag = flowSolver::HLLC;
	else if (name == "HLL")
		flowSolwerFlag = flowSolver::HLL;
	else if (name == "HLLCLM")
		flowSolwerFlag = flowSolver::HLLCLM;
	else if (name == "KT")
		flowSolwerFlag = flowSolver::KT;
	else if (name == "HLLC2p")
		flowSolwerFlag = flowSolver::HLLC2p;
	else if (name == "Richtmyer")
		flowSolwerFlag = flowSolver::Richtmyer;
	else
		throw exception("Unknown flow solver flag.",
				errors::initializationError);

	switch (flowSolwerFlag)
	{
	case flowSolver::HLL:
		return std::make_unique<HLLSolver>();
		break;
	case flowSolver::HLLCF:
		return std::make_unique<HLLCFSolver>();
		break;
	case flowSolver::HLLC:
		return std::make_unique<HLLCSolver>();
		break;
	case flowSolver::HLLCLM:
		return std::make_unique<HLLCLMSolver>();
		break;
	case flowSolver::HLLC2p:
		return std::make_unique<HLLC2pSolver>();
		break;
	case flowSolver::Richtmyer:
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
				errors::initializationError);

	turbulenceModel turbulenceModelFlag;
	if (sourceTypeString == "BHR")
		turbulenceModelFlag = turbulenceModel::BHRSource;
	else if (sourceTypeString == "zero")
		turbulenceModelFlag = turbulenceModel::zeroSource;
	else if (sourceTypeString == "decay")
		turbulenceModelFlag = turbulenceModel::decaySource;
	else if (sourceTypeString == "shear")
		turbulenceModelFlag = turbulenceModel::shearSource;
	else if (sourceTypeString == "arithmetic1")
		turbulenceModelFlag = turbulenceModel::arithmeticA1Source;
	else if (sourceTypeString == "arithmetic2")
		turbulenceModelFlag = turbulenceModel::arithmeticA2Source;
	else if (sourceTypeString == "kEpsA")
		turbulenceModelFlag = turbulenceModel::kEpsASource;
	else if (sourceTypeString == "arithmetic3")
		turbulenceModelFlag = turbulenceModel::arithmeticA3Source;
	else if (sourceTypeString == "BHRKL")
		turbulenceModelFlag = turbulenceModel::BHRKLSource;
	else
		throw exception("Unknown source generation flag.",
				errors::initializationError);

	std::unique_ptr<abstractTurbulenceGen> turbulenceSources;
	switch (turbulenceModelFlag)
	{
	case turbulenceModel::BHRSource:
		turbulenceSources = std::make_unique<BHRGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::decaySource:
		turbulenceSources = std::make_unique<decayGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::shearSource:
		turbulenceSources = std::make_unique<shearGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::arithmeticA1Source:
	case turbulenceModel::arithmeticA2Source:
	case turbulenceModel::arithmeticA3Source:
		turbulenceSources = std::make_unique<arithmeticAGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::kEpsASource:
		turbulenceSources = std::make_unique<kEpsAGen>(turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::BHRKLSource:
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
					errors::initializationError);

		turbulenceSources->turbPar->readTurbulentParameters(
				turbulentParametersFile);

		turbulentParametersFile.close();
	}

	return turbulenceSources;
}

std::unique_ptr<schemi::abstractChemicalKinetics> schemi::createChemicalKinetics(
		const homogeneousPhase<cubicCell> & phaseIn,
		const chemicalReactions chemReactFlag, const scalar minimalTimestep)
{
	switch (chemReactFlag)
	{
	case chemicalReactions::Cl2Dissociation:
		return std::make_unique<chemicalKineticsChlorumDissociation>(phaseIn,
				minimalTimestep);
		break;
	case chemicalReactions::Cl2H2Dissociation:
		return std::make_unique<chemicalKineticsChlorumHydrogeniumDissociation>(
				phaseIn, minimalTimestep);
		break;
	case chemicalReactions::H2Cl2Combustion:
		return std::make_unique<chemicalKineticsH2Cl2Combustion>(phaseIn,
				minimalTimestep);
		break;
	case chemicalReactions::NO2Disproportionation:
		return std::make_unique<chemicalKineticsNO2Disproportionation>(phaseIn,
				minimalTimestep);
		break;
	case chemicalReactions::noReaction:
	default:
		return std::make_unique<chemicalKineticsNoReaction>();
		break;
	}
}

std::unique_ptr<schemi::abstractTransportModel> schemi::createTransportModel(
		const std::vector<std::vector<std::string>> & matrixOfSubstancesConditions,
		const scalar constNu, const scalar constD, const scalar constKappa,
		const transportModel model)
{
	switch (model)
	{
	case transportModel::hardSpheres:
	{
		std::valarray<scalar> molDiams(matrixOfSubstancesConditions.size());

		for (std::size_t k = 0; k < molDiams.size(); ++k)
			molDiams[k] = std::stod(matrixOfSubstancesConditions[k][4]);

		return std::make_unique<hardSpheresTransportModel>(constNu, constD,
				constKappa, molDiams);
	}
		break;
	case transportModel::constant:
	default:
		return std::make_unique<abstractTransportModel>(constNu, constD,
				constKappa);
		break;
	}
}
