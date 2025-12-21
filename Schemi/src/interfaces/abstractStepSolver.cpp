/*
 * abstractStepSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractStepSolver.hpp"

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
		const boundaryConditionValue & bncCalc_in,
		const volumeField<scalar> & minimalLengthScale_in,
		const timestep & sourceTimeFlag_in,
		const bool & molMassDiffusionFlag_in,
		const chemicalKinetics::abstractChemicalKinetics & chemKin_in,
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
				enthalpyFlowFlag_in), linearFlag(linearFlag_in), bncCalc(
				bncCalc_in), minimalLengthScale(minimalLengthScale_in), sourceTimeFlag(
				sourceTimeFlag_in), molMassDiffusionFlag(
				molMassDiffusionFlag_in), chemKin(chemKin_in), nonLinearityIteratonsFlag(
				nonLinearityIteratonsFlag_in)
{
}

schemi::abstractStepSolver::~abstractStepSolver() noexcept
{
}
