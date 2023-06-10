/*
 * HLLC2pSolver.hpp
 *
 *  Created on: 2021/12/18
 *      Author: Maxim Boldyrev
 *
 *      HLLC Riemann with two star pressures solver class.
 */

#ifndef HLLC2PSOLVER_HPP_
#define HLLC2PSOLVER_HPP_

#include "abstractFlowSolver.hpp"
#include "pressureStarClass.hpp"

namespace schemi
{
class HLLC2pSolver: public abstractFlowSolver, private pressureStarClass
{
public:
	std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
			const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
					override;
};
}  // namespace schemi

#endif /* HLLC2PSOLVER_HPP_ */
