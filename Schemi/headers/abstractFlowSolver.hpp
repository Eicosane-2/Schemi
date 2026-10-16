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
#include "homogeneousPhase.hpp"
#include "MPIHandler.hpp"
#include "starFields.hpp"

namespace schemi
{
class abstractFlowSolver
{
protected:
	std::array<scalar, 2> DavisSignalSpeeds(const scalar uL, const scalar aL,
			const scalar uR, const scalar aR) const noexcept;
	std::array<scalar, 2> EinfeldtSignalSpeeds(const scalar rhoL,
			const scalar uL, const scalar aL, const scalar rhoR,
			const scalar uR, const scalar aR) const noexcept;
public:
	virtual ~abstractFlowSolver() noexcept =0;

	static std::unique_ptr<abstractFlowSolver> createFlowSolver(
			const std::string & name, const MPIHandler & par);

	virtual std::tuple<conservativeFlows, starFields> calculateFlows(
			const homogeneousPhase<quadraticSurface>& /*surfaceOwnerSide*/,
			const homogeneousPhase<quadraticSurface>& /*surfaceNeighbourSide*/) const =0;
};
}  // namespace schemi

#endif /* ABSTRACTFLOWSOLVER_HPP_ */
