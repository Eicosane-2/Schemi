/*
 * abstractStepSolver.hpp
 *
 *  Created on: 2022/12/11
 *      Author: Maxim Boldyrev
 */

#ifndef ABSTRACTSTEPSOLVER_HPP_
#define ABSTRACTSTEPSOLVER_HPP_

#include <vector>

#include "boundaryConditionTypesEnum.hpp"
#include "enthalpyFlowEnum.hpp"
#include "timestepEnum.hpp"
#include "abstractChemicalKinetics.hpp"
#include "abstractFlowSolver.hpp"
#include "abstractLimiter.hpp"
#include "abstractMatrixSolver.hpp"
#include "boundaryConditionValue.hpp"
#include "cubicCell.hpp"
#include "MPIHandler.hpp"
#include "scalar.hpp"
#include "homogeneousPhase.hpp"
#include "volumeField.hpp"

namespace schemi
{
class abstractStepSolver
{
protected:
	homogeneousPhase<cubicCell> & gasPhase;
	const abstractLimiter & limiter;
	const abstractFlowSolver & fsolver;
	const bool & gravitationFlag;
	const vector & g;
	const boundaryConditionValue & boundaryConditionValueCalc;
	scalar & timeForTVD;
	scalar & timeForHancock;
	scalar & timeForFlowCalculation;
	scalar & timeForTimeIntegration;
	const MPIHandler & parallelism;
	const bool & diffusionFlag;
	const abstractMatrixSolver & msolver;
	const abstractMatrixSolver & msolverEnthFl;
	const std::pair<scalar, scalar> timestepCoeffs;
	scalar & timeForDiffusion;
	const std::vector<boundaryConditionType> & commonConditions;
	const enthalpyFlow & enthalpyFlowFlag;
	const bool & linearFlag;
	const boundaryConditionValue & bncCalc;
	const volumeField<scalar> & minimalLengthScale;
	const timestep & sourceTimeFlag;
	const bool & molMassDiffusionFlag;
	const abstractChemicalKinetics & chemKin;
public:
	abstractStepSolver(homogeneousPhase<cubicCell> & gasPhase_in,
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
			const abstractChemicalKinetics & chemKin_in) noexcept;

	abstractStepSolver(const abstractStepSolver&) = delete;
	auto& operator=(const abstractStepSolver&) = delete;

	virtual ~abstractStepSolver() noexcept =0;

	virtual void calculateStep() =0;
};
}  // namespace schemi

#endif /* ABSTRACTSTEPSOLVER_HPP_ */
