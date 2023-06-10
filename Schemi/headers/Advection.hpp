/*
 * Advection.hpp
 *
 *  Created on: 2020/04/29
 *      Author: Maxim Boldyrev
 *
 *      Function for solving advection system of equations.
 */

#ifndef ADVECTION_HPP_
#define ADVECTION_HPP_

#include "homogeneousPhase.hpp"
#include "abstractFlowSolver.hpp"
#include "abstractLimiter.hpp"
#include "boundaryConditionValue.hpp"
#include "cubicCell.hpp"
#include "MPIHandler.hpp"
#include "scalar.hpp"
#include "starFields.hpp"

namespace schemi
{
/*Advection stage.*/
starFields Advection(homogeneousPhase<cubicCell> & gasPhase,
		const abstractLimiter & limiter, const abstractFlowSolver & fsolver,
		std::pair<bool, vector> gravitation,
		const boundaryConditionValue & boundaryConditionValueCalc,
		scalar & timeForTVD, scalar & timeForHancock,
		scalar & timeForFlowCalculation, scalar & timeForTimeIntegration,
		const MPIHandler & parallelism);
}  // namespace schemi

#endif /* ADVECTION_HPP_ */
