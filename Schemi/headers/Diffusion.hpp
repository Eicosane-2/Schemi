/*
 * Diffusion.hpp
 *
 *  Created on: 2020/04/29
 *      Author: Maxim Boldyrev
 *
 *      Function for solving diffusion and turbulence sources system of equations.
 */

#ifndef DIFFUSION_HPP_
#define DIFFUSION_HPP_

#include "enthalpyFlowEnum.hpp"
#include "timestepEnum.hpp"
#include "abstractMatrixSolver.hpp"
#include "cubicCell.hpp"
#include "MPIHandler.hpp"
#include "starFields.hpp"
#include "homogeneousPhase.hpp"

namespace schemi
{
/*Diffusion stage.*/
void Diffusion(homogeneousPhase<cubicCell> & gasPhase,
		const abstractMatrixSolver & msolver,
		const abstractMatrixSolver & msolverForEnthalpy,
		const std::pair<scalar, scalar> & timestepCoeffs,
		scalar & timeForDiffusion,
		const std::vector<boundaryConditionType> & commonConditions,
		const starFields & star, const enthalpyFlowEnum enthalpySolverFlag,
		const bool linearRec, const boundaryConditionValue & bncCalc,
		const volumeField<scalar> & minimalLengthScale,
		const MPIHandler & parallelism, const timestepEnum sourceTimeFlag,
		const bool molMassDiffusionFlag);
}  // namespace schemi

#endif /* DIFFUSION_HPP_ */
