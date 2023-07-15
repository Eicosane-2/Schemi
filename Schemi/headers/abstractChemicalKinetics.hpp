/*
 * abstractChemicalKinetics.hpp
 *
 *  Created on: 2023/05/08
 *      Author: Maxim Boldyrev
 */

#ifndef ABSTRACTCHEMICALKINETICS_HPP_
#define ABSTRACTCHEMICALKINETICS_HPP_

#include <cstddef>
#include <vector>
#include <utility>

#include "cubicCell.hpp"
#include "scalar.hpp"
#include "homogeneousPhase.hpp"

namespace schemi
{
class abstractChemicalKinetics
{
protected:
	typedef std::vector<std::pair<scalar, std::size_t>> triangleList;

	constexpr static scalar convergenceTolerance { 1E-12 };
	std::size_t maxIterationNumber { 0 };

	enum class iterativeSolver
	{
		noSolver, GaussSeidel, ConjugateGradient, JacobiConjugateGradient
	};
public:
	const bool chemicalReaction;

	abstractChemicalKinetics(const bool flag) noexcept;

	virtual ~abstractChemicalKinetics() noexcept =0;

	virtual void solveChemicalKinetics(
			homogeneousPhase<cubicCell>&) const noexcept =0;
};
}  // namespace schemi

#endif /* ABSTRACTCHEMICALKINETICS_HPP_ */
