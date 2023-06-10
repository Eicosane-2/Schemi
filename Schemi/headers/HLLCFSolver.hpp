/*
 * HLLCFSolver.hpp
 *
 *  Created on: 2019/12/23
 *      Author: Maxim Boldyrev
 *
 *      HLLC-Focus Riemann solver class.
 */

#ifndef HLLCFSOLVER_HPP_
#define HLLCFSOLVER_HPP_

#include "abstractFlowSolver.hpp"
#include "pressureStarClass.hpp"

namespace schemi
{
class HLLCFSolver: public abstractFlowSolver, private pressureStarClass
{
public:
	std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
			const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
					override;
};
}  // namespace schemi

#endif /* HLLCFSOLVER_HPP_ */
