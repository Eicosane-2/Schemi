/*
 * HLLCLMSolver.hpp
 *
 *  Created on: 2021/02/08
 *      Author: Maxim Boldyrev
 *
 *      HLLC-Low Mach Riemann solver class.
 */

#ifndef HLLCLMSOLVER_HPP_
#define HLLCLMSOLVER_HPP_

#include "abstractFlowSolver.hpp"
#include "pressureStarClass.hpp"

namespace schemi
{
class HLLCLMSolver: public abstractFlowSolver, private pressureStarClass
{
	constexpr static scalar minMachNumber = 0.1;

public:
	std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
			const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
					override;
};
}  // namespace schemi

#endif /* HLLCLMSOLVER_HPP_ */
