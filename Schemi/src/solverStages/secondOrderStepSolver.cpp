/*
 * secondOrderStepSolver.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "secondOrderStepSolver.hpp"

#include "Advection.hpp"
#include "Diffusion.hpp"
#include "linearInterpolate.hpp"

schemi::secondOrderStepSolver::secondOrderStepSolver(
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
		const chemicalKinetics::abstractChemicalKinetics & chemKin_in) noexcept :
		abstractStepSolver(gasPhase_in, limiter_in, fsolver_in,
				gravitationFlag_in, g_in, boundaryConditionValueCalc_in,
				timeForTVD_in, timeForHancock_in, timeForFlowCalculation_in,
				timeForTimeIntegration_in, parallelism_in, diffusionFlag_in,
				msolver_in, msolverEnthFl_in, timestepCoeffs_in,
				timeForDiffusion_in, commonConditions_in, enthalpyFlowFlag_in,
				linearFlag_in, bncCalc_in, minimalLengthScale_in,
				sourceTimeFlag_in, molMassDiffusionFlag_in, chemKin_in)
{
}

void schemi::secondOrderStepSolver::calculateStep()
{
	auto star = Advection(gasPhase, limiter, fsolver, { gravitationFlag, g },
			boundaryConditionValueCalc, timeForTVD, timeForHancock,
			timeForFlowCalculation, timeForTimeIntegration, parallelism);

	parallelism.correctBoundaryValues(gasPhase);

	if (diffusionFlag)
	{
		if (gravitationFlag) /*Gravitational acceleration invalidates star values of velocity.*/
			star.v = linearInterpolate(gasPhase.velocity,
					boundaryConditionValueCalc);

		Diffusion(gasPhase, msolver, msolverEnthFl, timestepCoeffs,
				timeForDiffusion, commonConditions, star, enthalpyFlowFlag,
				linearFlag, boundaryConditionValueCalc, minimalLengthScale,
				parallelism, sourceTimeFlag, molMassDiffusionFlag);

		parallelism.correctBoundaryValues(gasPhase);
	}

	if (chemKin.chemicalReaction)
	{
		chemKin.solveChemicalKinetics(gasPhase);

		parallelism.correctBoundaryValues(gasPhase);
	}
}
