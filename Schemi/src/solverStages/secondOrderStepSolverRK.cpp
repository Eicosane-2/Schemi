/*
 * secondOrderStepSolverRK.cpp
 *
 *  Created on: 2024/11/20
 *      Author: Maxim Boldyrev
 */

#include "secondOrderStepSolverRK.hpp"

#include "Advection2dOrder.hpp"
#include "Diffusion.hpp"
#include "linearInterpolate.hpp"
#include "gradient.hpp"
#include "divergence.hpp"

schemi::secondOrderStepSolverRK::secondOrderStepSolverRK(
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
		chemicalKinetics::abstractChemicalKinetics & chemKin_in,
		const bool & nonLinearityIteratonsFlag_in) noexcept :
		abstractStepSolver(gasPhase_in, limiter_in, fsolver_in,
				gravitationFlag_in, g_in, boundaryConditionValueCalc_in,
				timeForTVD_in, timeForHancock_in, timeForFlowCalculation_in,
				timeForTimeIntegration_in, parallelism_in, diffusionFlag_in,
				msolver_in, msolverEnthFl_in, timestepCoeffs_in,
				timeForDiffusion_in, commonConditions_in, enthalpyFlowFlag_in,
				linearFlag_in, bncCalc_in, minimalLengthScale_in,
				sourceTimeFlag_in, molMassDiffusionFlag_in, chemKin_in,
				nonLinearityIteratonsFlag_in)
{
}

void schemi::secondOrderStepSolverRK::calculateStep()
{
	const bunchOfFields<cubicCell> Un = gasPhase;

	auto star = Advection2dOrder(gasPhase, limiter, fsolver, { gravitationFlag,
			g }, boundaryConditionValueCalc, timeForTVD, timeForHancock,
			timeForFlowCalculation, timeForTimeIntegration, parallelism);

	/* gasPhase: Un1 = Un - dt*divF(Un) */

	parallelism.correctBoundaryValues(gasPhase);

	const auto star1 = Advection2dOrder(gasPhase, limiter, fsolver, {
			gravitationFlag, g }, boundaryConditionValueCalc, timeForTVD,
			timeForHancock, timeForFlowCalculation, timeForTimeIntegration,
			parallelism);

	/* gasPhase: Un2 = Un1 - dt*divF(Un1) */

	bunchOfFields<cubicCell> Un2 = gasPhase;
	Un2.average(Un, *gasPhase.phaseThermodynamics, 0.5);
	gasPhase.copyFrom(Un2, *gasPhase.phaseThermodynamics);

	parallelism.correctBoundaryValues(gasPhase);

	star.c.v[0].val() = 0;
	for (std::size_t k = 1; k < star.c.v.size(); ++k)
	{
		star.c.v[k] = (star.c.v[k] + star1.c.v[k]) / 2;

		star.c.v[0] += star.c.v[k];
	}

	star.rho = (star.rho + star1.rho) / 2;
	star.v = (star.v + star1.v) / 2;
	star.p = (star.p + star1.p) / 2;
	star.a = (star.a + star1.a) / 2;
	star.b = (star.b + star1.b) / 2;

	if (gasPhase.turbulence->isInitialisationModelUsed())
	{
		const volumeField<vector> gradRho_c(grad(star.rho));
		const surfaceField<vector> gradRho_s(
				surfGrad(gasPhase.density[0], bncCalc));
		const volumeField<vector> gradP(grad(star.p));
		const volumeField<scalar> divU(divergence(star.v));
		const volumeField<tensor> gradU(grad(star.v));

		gasPhase.turbulence->particlesTimeIntegration(gradRho_c, gradRho_s,
				gasPhase.velocity, star.v, g * gravitationFlag,
				gasPhase.concentration, gasPhase.density, bncCalc,
				gasPhase.phaseThermodynamics->Mv(),
				gasPhase.pressure.meshRef().timestep(), gradP, divU, gradU);
	}

	if (diffusionFlag)
	{
		Diffusion(gasPhase, msolver, msolverEnthFl, timestepCoeffs,
				timeForDiffusion, commonConditions, star, enthalpyFlowFlag,
				linearFlag, boundaryConditionValueCalc, minimalLengthScale,
				parallelism, sourceTimeFlag, molMassDiffusionFlag,
				nonLinearityIteratonsFlag);

		parallelism.correctBoundaryValues(gasPhase);
	}

	if (chemKin.chemicalReaction)
	{
		chemKin.solveChemicalKinetics(gasPhase);

		parallelism.correctBoundaryValues(gasPhase);
	}
}
