/*
 * HLLCKTSolver.hpp
 *
 *  Created on: 2026/09/11
 *      Author: Maxim Boldyrev
 */

#ifndef HLLCKTSOLVER_HPP_
#define HLLCKTSOLVER_HPP_

#include "abstractFlowSolver.hpp"
#include "pressureStarClass.hpp"

namespace schemi
{
class HLLCKTSolver: public abstractFlowSolver, private pressureStarClass
{
public:
	std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
			const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
					override;
};
}  // namespace schemi

#endif /* HLLCKTSOLVER_HPP_ */
