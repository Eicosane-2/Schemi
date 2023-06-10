/*
 * RichtmyerSolver.hpp
 *
 *  Created on: 2023/01/08
 *      Author: Maxim Boldyrev
 */

#ifndef RICHTMYERSOLVER_HPP_
#define RICHTMYERSOLVER_HPP_

#include "abstractFlowSolver.hpp"
#include "pressureStarClass.hpp"

namespace schemi
{
class RichtmyerSolver: public abstractFlowSolver, private pressureStarClass
{
public:
	std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
			const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
					override;
};
}  // namespace schemi

#endif /* RICHTMYERSOLVER_HPP_ */
