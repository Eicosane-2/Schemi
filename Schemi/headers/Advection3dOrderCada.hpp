/*
 * Advection3dOrderCada.hpp
 *
 *  Created on: 2025/06/05
 *      Author: Maxim Boldyrev
 */

#ifndef ADVECTION3DORDERCADA_HPP_
#define ADVECTION3DORDERCADA_HPP_

#include "abstractFlowSolver.hpp"
#include "abstractLimiter.hpp"
#include "boundaryConditionValue.hpp"
#include "MPIHandler.hpp"
#include "starFields.hpp"
#include "homogeneousPhase.hpp"

namespace schemi
{
/*Advection stage.*/
starFields Advection3dOrderCada(homogeneousPhase<cubicCell> & gasPhase,
		const abstractLimiter & limiter, const abstractFlowSolver & fsolver,
		std::pair<bool, vector> gravitation,
		const boundaryConditionValue & boundaryConditionValueCalc,
		scalar & timeForTVD, scalar & timeForHancock,
		scalar & timeForFlowCalculation, scalar & timeForTimeIntegration,
		const MPIHandler & parallelism,
		const volumeField<scalar> & minimalLengthScale);
}  // namespace schemi

#endif /* ADVECTION3DORDERCADA_HPP_ */
