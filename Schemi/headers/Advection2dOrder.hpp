/*
 * Advection2dOrder.hpp
 *
 *  Created on: 2024/11/20
 *      Author: Maxim Boldyrev
 */

#ifndef ADVECTION2DORDER_HPP_
#define ADVECTION2DORDER_HPP_

#include "abstractFlowSolver.hpp"
#include "abstractLimiter.hpp"
#include "boundaryConditionValue.hpp"
#include "MPIHandler.hpp"
#include "starFields.hpp"
#include "homogeneousPhase.hpp"

namespace schemi
{
/*Advection stage.*/
starFields Advection2dOrder(homogeneousPhase<cubicCell> & gasPhase,
		const abstractLimiter & limiter, const abstractFlowSolver & fsolver,
		std::pair<bool, vector> gravitation,
		const boundaryConditionValue & boundaryConditionValueCalc,
		scalar & timeForTVD, scalar & timeForHancock,
		scalar & timeForFlowCalculation, scalar & timeForTimeIntegration,
		const MPIHandler & parallelism);
}  // namespace schemi

#endif /* ADVECTION2DORDER_HPP_ */
