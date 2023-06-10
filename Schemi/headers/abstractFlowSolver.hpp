/*
 * abstractFlowSolver.hpp
 *
 *  Created on: 2020/03/31
 *      Author: Maxim Boldyrev
 *
 *      Interface class for Riemann solvers.
 */

#ifndef ABSTRACTFLOWSOLVER_HPP_
#define ABSTRACTFLOWSOLVER_HPP_

#include <tuple>

#include "conservativeFlows.hpp"
#include "quadraticSurface.hpp"
#include "starFields.hpp"
#include "homogeneousPhase.hpp"

namespace schemi
{
class abstractFlowSolver
{
public:
	virtual ~abstractFlowSolver() noexcept =0;

	virtual std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface>& /*surfaceOwnerSide*/,
			const homogeneousPhase<quadraticSurface>& /*surfaceNeighbourSide*/) const =0;
};
}  // namespace schemi

#endif /* ABSTRACTFLOWSOLVER_HPP_ */
