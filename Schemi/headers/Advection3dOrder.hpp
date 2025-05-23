/*
 * Advection3dOrder.hpp
 *
 *  Created on: 2022/12/03
 *      Author: Maxim Boldyrev
 *
 *      Function for solving advection system of equations with 3d order spatial approximation.
 */

#ifndef ADVECTION3DORDER_HPP_
#define ADVECTION3DORDER_HPP_

#include "abstractFlowSolver.hpp"
#include "abstractLimiter.hpp"
#include "boundaryConditionValue.hpp"
#include "MPIHandler.hpp"
#include "starFields.hpp"
#include "homogeneousPhase.hpp"

namespace schemi
{
/*Advection stage.*/
starFields Advection3dOrder(homogeneousPhase<cubicCell> & gasPhase,
		const abstractLimiter & limiter, const abstractFlowSolver & fsolver,
		std::pair<bool, vector> gravitation,
		const boundaryConditionValue & boundaryConditionValueCalc,
		scalar & timeForTVD, scalar & timeForHancock,
		scalar & timeForFlowCalculation, scalar & timeForTimeIntegration,
		const MPIHandler & parallelism,
		const volumeField<scalar> & minimalLengthScale);
}  // namespace schemi

#endif /* ADVECTION3DORDER_HPP_ */
