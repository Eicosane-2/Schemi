/*
 * KTSolver.hpp
 *
 *  Created on: 2021/03/12
 *      Author: Maxim Boldyrev
 *
 *      Kurganov-Tadmor solver class.
 */

#ifndef KTSOLVER_HPP_
#define KTSOLVER_HPP_

#include <tuple>

#include "abstractFlowSolver.hpp"
#include "pressureStarClass.hpp"

namespace schemi
{
class KTSolver: public abstractFlowSolver, private pressureStarClass
{
public:
	std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
			const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
					override;
};
}  // namespace schemi

#endif /* KTSOLVER_HPP_ */
