/*
 * thirdOrderStepSolverCada.hpp
 *
 *  Created on: 2025/06/05
 *      Author: Maxim Boldyrev
 */

#ifndef THIRDORDERSTEPSOLVERCADA_HPP_
#define THIRDORDERSTEPSOLVERCADA_HPP_

#include "abstractStepSolver.hpp"

namespace schemi
{
class thirdOrderStepSolverCada: public abstractStepSolver
{
public:
	thirdOrderStepSolverCada(homogeneousPhase<cubicCell> & gasPhase_in,
			const abstractLimiter & limiter_in,
			const abstractFlowSolver & fsolver_in,
			const bool & gravitationFlag_in, const vector & g_in,
			const boundaryConditionValue & boundaryConditionValueCalc_in,
			scalar & timeForTVD_in, scalar & timeForHancock_in,
			scalar & timeForFlowCalculation_in,
			scalar & timeForTimeIntegration_in,
			const MPIHandler & parallelism_in, const bool & diffusionFlag_in,
			const abstractMatrixSolver & msolver_in,
			const abstractMatrixSolver & msolverEnthFl_in,
			const std::pair<scalar, scalar> & timestepCoeffs_in,
			scalar & timeForDiffusion_in,
			const std::vector<boundaryConditionType> & commonConditions_in,
			const enthalpyFlow & enthalpyFlowFlag_in,
			const bool & linearFlag_in,
			const boundaryConditionValue & bncCalc_in,
			const volumeField<scalar> & minimalLengthScale_in,
			const timestep & sourceTimeFlag_in,
			const bool & molMassDiffusionFlag_in,
			chemicalKinetics::abstractChemicalKinetics & chemKin_in,
			const bool & nonLinearityIteratonsFlag_in) noexcept;

	thirdOrderStepSolverCada(const thirdOrderStepSolverCada&) = delete;
	auto& operator=(const thirdOrderStepSolverCada&) = delete;

	void calculateStep() override;
};
}  // namespace schemi

#endif /* THIRDORDERSTEPSOLVERCADA_HPP_ */
