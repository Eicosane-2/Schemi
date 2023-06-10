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

#include "homogeneousPhase.hpp"
#include "cubicCell.hpp"
#include "scalar.hpp"

namespace schemi
{
class abstractChemicalKinetics
{
protected:
	constexpr static scalar convergenceTolerance { 1E-12 };
	const std::size_t maxIterationNumber;
	typedef std::vector<std::pair<scalar, std::size_t>> triangleList;
public:
	const bool chemicalReaction;

	abstractChemicalKinetics(const bool flag, const std::size_t number) noexcept;

	virtual ~abstractChemicalKinetics() noexcept =0;

	virtual void solveChemicalKinetics(
			homogeneousPhase<cubicCell>&) const noexcept =0;
};
}  // namespace schemi

#endif /* ABSTRACTCHEMICALKINETICS_HPP_ */
