/*
 * thirdOrderStepSolver.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "thirdOrderStepSolver.hpp"

#include "Advection3dOrder.hpp"
#include "Diffusion.hpp"
#include "linearInterpolate.hpp"

schemi::thirdOrderStepSolver::thirdOrderStepSolver(
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
		const abstractChemicalKinetics & chemKin_in) noexcept :
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

void schemi::thirdOrderStepSolver::calculateStep()
{
	const bunchOfFields<cubicCell> Un = gasPhase;

	auto star = Advection3dOrder(gasPhase, limiter, fsolver, { gravitationFlag,
			g }, boundaryConditionValueCalc, timeForTVD, timeForHancock,
			timeForFlowCalculation, timeForTimeIntegration, parallelism);

	/* gasPhase: Un1 */

	parallelism.correctBoundaryValues(gasPhase);

	Advection3dOrder(gasPhase, limiter, fsolver, { gravitationFlag, g },
			boundaryConditionValueCalc, timeForTVD, timeForHancock,
			timeForFlowCalculation, timeForTimeIntegration, parallelism);

	/* gasPhase: Un1 - dt*divF(Un1) */

	bunchOfFields<cubicCell> Un2 = gasPhase;
	Un2.average(Un, *gasPhase.phaseThermodynamics, 0.25);
	gasPhase.copyFrom(Un2, *gasPhase.phaseThermodynamics);

	parallelism.correctBoundaryValues(gasPhase);

	Advection3dOrder(gasPhase, limiter, fsolver, { gravitationFlag, g },
			boundaryConditionValueCalc, timeForTVD, timeForHancock,
			timeForFlowCalculation, timeForTimeIntegration, parallelism);

	/* gasPhase: Un2 - dt*divF(Un2) */

	bunchOfFields<cubicCell> Un3 = gasPhase;
	Un3.average(Un, *gasPhase.phaseThermodynamics, 2. / 3.);
	gasPhase.copyFrom(Un3, *gasPhase.phaseThermodynamics);

	parallelism.correctBoundaryValues(gasPhase);

	if (diffusionFlag)
	{
		star.c.v[0].ref_r() = 0;
		for (std::size_t k = 1; k < star.c.v.size(); ++k)
		{
			star.c.v[k] = linearInterpolate(gasPhase.concentration.v[k],
					boundaryConditionValueCalc);

			star.c.v[0].ref_r() += star.c.v[k].ref();
		}

		star.rho = linearInterpolate(gasPhase.density[0],
				boundaryConditionValueCalc);
		star.v = linearInterpolate(gasPhase.velocity,
				boundaryConditionValueCalc);
		star.p = linearInterpolate(gasPhase.pressure,
				boundaryConditionValueCalc);
		star.a = linearInterpolate(gasPhase.aTurb, boundaryConditionValueCalc);
		star.b = linearInterpolate(gasPhase.bTurb, boundaryConditionValueCalc);

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
