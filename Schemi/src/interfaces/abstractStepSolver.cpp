/*
 * abstractStepSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractStepSolver.hpp"

#include "typeOfSolverEnum.hpp"
#include "secondOrderStepSolver.hpp"
#include "secondOrderStepSolverRK.hpp"
#include "thirdOrderStepSolver.hpp"
#include "thirdOrderStepSolverCada.hpp"

schemi::abstractStepSolver::abstractStepSolver(
		homogeneousPhase<cubicCell> & gasPhase_in,
		const abstractLimiter & limiter_in,
		const abstractFlowSolver & fsolver_in, const bool & gravitationFlag_in,
		const vector & g_in,
		const boundaryConditionValue & boundaryConditionValueCalc_in,
		scalar & timeForTVD_in, scalar & timeForHancock_in,
		scalar & timeForFlowCalculation_in, scalar & timeForTimeIntegration_in,
		const MPIHandler & parallelism_in, const bool & diffusionFlag_in,
		const abstractMatrixSolver & msolver_in,
		const abstractMatrixSolver & msolverEnthFl_in,
		const std::pair<scalar, scalar> & timestepCoeffs_in,
		scalar & timeForDiffusion_in,
		const std::vector<boundaryConditionType> & commonConditions_in,
		const enthalpyFlow & enthalpyFlowFlag_in, const bool & linearFlag_in,
		const volumeField<scalar> & minimalLengthScale_in,
		const timestep & sourceTimeFlag_in,
		const bool & molMassDiffusionFlag_in,
		chemicalKinetics::abstractChemicalKinetics & chemKin_in,
		const bool & nonLinearityIteratonsFlag_in) noexcept :
		gasPhase(gasPhase_in), limiter(limiter_in), fsolver(fsolver_in), gravitationFlag(
				gravitationFlag_in), g(g_in), boundaryConditionValueCalc(
				boundaryConditionValueCalc_in), timeForTVD(timeForTVD_in), timeForHancock(
				timeForHancock_in), timeForFlowCalculation(
				timeForFlowCalculation_in), timeForTimeIntegration(
				timeForTimeIntegration_in), parallelism(parallelism_in), diffusionFlag(
				diffusionFlag_in), msolver(msolver_in), msolverEnthFl(
				msolverEnthFl_in), timestepCoeffs(timestepCoeffs_in), timeForDiffusion(
				timeForDiffusion_in), commonConditions(commonConditions_in), enthalpyFlowFlag(
				enthalpyFlowFlag_in), linearFlag(linearFlag_in), minimalLengthScale(
				minimalLengthScale_in), sourceTimeFlag(sourceTimeFlag_in), molMassDiffusionFlag(
				molMassDiffusionFlag_in), chemKin(chemKin_in), nonLinearityIteratonsFlag(
				nonLinearityIteratonsFlag_in)
{
}

schemi::abstractStepSolver::~abstractStepSolver() noexcept
{
}

std::unique_ptr<schemi::abstractStepSolver> schemi::abstractStepSolver::createStepSolver(
		const std::string & thirdOrderString,
		homogeneousPhase<cubicCell> & gasPhase_in,
		const abstractLimiter & limiter_in,
		const abstractFlowSolver & fsolver_in, const bool & gravitationFlag_in,
		const vector & g_in,
		const boundaryConditionValue & boundaryConditionValueCalc_in,
		scalar & timeForTVD_in, scalar & timeForHancock_in,
		scalar & timeForFlowCalculation_in, scalar & timeForTimeIntegration_in,
		const MPIHandler & parallelism_in, const bool & diffusionFlag_in,
		const abstractMatrixSolver & msolver_in,
		const abstractMatrixSolver & msolverEnthFl_in,
		const std::pair<scalar, scalar> & timestepCoeffs_in,
		scalar & timeForDiffusion_in,
		const std::vector<boundaryConditionType> & commonConditions_in,
		const enthalpyFlow & enthalpyFlowFlag_in, const bool & linearFlag_in,
		const volumeField<scalar> & minimalLengthScale_in,
		const timestep & sourceTimeFlag_in,
		const bool & molMassDiffusionFlag_in,
		chemicalKinetics::abstractChemicalKinetics & chemKin_in,
		const bool & nonLinearityIteratonsFlag_in)
{
	typeOfSolverEnum order;

	std::map<std::string, typeOfSolverEnum> orderType;
	orderType.insert( { "ThirdOrder", typeOfSolverEnum::ThirdOrder });
	orderType.insert( { "ThirdOrderCada", typeOfSolverEnum::ThirdOrderCada });
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

	switch (order)
	{
	case typeOfSolverEnum::SecondOrder:
		return std::make_unique<secondOrderStepSolver>(gasPhase_in, limiter_in,
				fsolver_in, gravitationFlag_in, g_in,
				boundaryConditionValueCalc_in, timeForTVD_in, timeForHancock_in,
				timeForFlowCalculation_in, timeForTimeIntegration_in,
				parallelism_in, diffusionFlag_in, msolver_in, msolverEnthFl_in,
				timestepCoeffs_in, timeForDiffusion_in, commonConditions_in,
				enthalpyFlowFlag_in, linearFlag_in, minimalLengthScale_in,
				sourceTimeFlag_in, molMassDiffusionFlag_in, chemKin_in,
				nonLinearityIteratonsFlag_in);
		break;
	case typeOfSolverEnum::SecondOrderRK:
		return std::make_unique<secondOrderStepSolverRK>(gasPhase_in,
				limiter_in, fsolver_in, gravitationFlag_in, g_in,
				boundaryConditionValueCalc_in, timeForTVD_in, timeForHancock_in,
				timeForFlowCalculation_in, timeForTimeIntegration_in,
				parallelism_in, diffusionFlag_in, msolver_in, msolverEnthFl_in,
				timestepCoeffs_in, timeForDiffusion_in, commonConditions_in,
				enthalpyFlowFlag_in, linearFlag_in, minimalLengthScale_in,
				sourceTimeFlag_in, molMassDiffusionFlag_in, chemKin_in,
				nonLinearityIteratonsFlag_in);
		break;
	case typeOfSolverEnum::ThirdOrder:
		return std::make_unique<thirdOrderStepSolver>(gasPhase_in, limiter_in,
				fsolver_in, gravitationFlag_in, g_in,
				boundaryConditionValueCalc_in, timeForTVD_in, timeForHancock_in,
				timeForFlowCalculation_in, timeForTimeIntegration_in,
				parallelism_in, diffusionFlag_in, msolver_in, msolverEnthFl_in,
				timestepCoeffs_in, timeForDiffusion_in, commonConditions_in,
				enthalpyFlowFlag_in, linearFlag_in, minimalLengthScale_in,
				sourceTimeFlag_in, molMassDiffusionFlag_in, chemKin_in,
				nonLinearityIteratonsFlag_in);
		break;
	case typeOfSolverEnum::ThirdOrderCada:
		return std::make_unique<thirdOrderStepSolverCada>(gasPhase_in,
				limiter_in, fsolver_in, gravitationFlag_in, g_in,
				boundaryConditionValueCalc_in, timeForTVD_in, timeForHancock_in,
				timeForFlowCalculation_in, timeForTimeIntegration_in,
				parallelism_in, diffusionFlag_in, msolver_in, msolverEnthFl_in,
				timestepCoeffs_in, timeForDiffusion_in, commonConditions_in,
				enthalpyFlowFlag_in, linearFlag_in, minimalLengthScale_in,
				sourceTimeFlag_in, molMassDiffusionFlag_in, chemKin_in,
				nonLinearityIteratonsFlag_in);
		break;
	[[unlikely]] default:
		throw exception("Unknown type of approximation order.",
				errors::initialisationError);
		break;
	}
}
