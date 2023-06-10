/*
 * HLLCSolver.hpp
 *
 *  Created on: 2020/04/01
 *      Author: Maxim Boldyrev
 *
 *      Classic HLLC Riemann solver class.
 */

#ifndef HLLCSOLVER_HPP_
#define HLLCSOLVER_HPP_

#include "abstractFlowSolver.hpp"
#include "pressureStarClass.hpp"

namespace schemi
{
class HLLCSolver: public abstractFlowSolver, private pressureStarClass
{
public:
	std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
			const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
					override;
};
}  // namespace schemi

#endif /* HLLCSOLVER_HPP_ */
