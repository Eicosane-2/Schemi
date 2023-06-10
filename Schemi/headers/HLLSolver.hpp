/*
 * HLLSolver.hpp
 *
 *  Created on: 2020/03/31
 *      Author: Maxim Boldyrev
 *
 *      HLL Riemann solver class.
 */

#ifndef HLLSOLVER_HPP_
#define HLLSOLVER_HPP_

#include "abstractFlowSolver.hpp"
#include "pressureStarClass.hpp"

namespace schemi
{
class HLLSolver: public abstractFlowSolver, private pressureStarClass
{
public:
	std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
			const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
					override;
};
}  // namespace schemi

#endif /* HLLSOLVER_HPP_ */
