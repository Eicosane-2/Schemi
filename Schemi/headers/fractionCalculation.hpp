/*
 * fractionCalculation.hpp
 *
 *  Created on: 2020/08/15
 *      Author: Maxim Boldyrev
 *
 *      Abstract class for mass and molar fraction calculation.
 */

#ifndef FRACTIONCALCULATION_HPP_
#define FRACTIONCALCULATION_HPP_

#include <valarray>
#include <vector>

#include "scalar.hpp"

namespace schemi
{
class fractionCalculation
{
protected:
	std::valarray<std::valarray<scalar>> calcMolarFrac(
			const std::vector<const std::valarray<scalar>*> & concentrations) const noexcept;

	std::valarray<std::valarray<scalar>> rearrangeMolFrac(
			const std::valarray<std::valarray<scalar>> & in) const noexcept;

	std::valarray<scalar> calcMolarFrac(
			const std::valarray<scalar> & concentrations) const noexcept;
public:
	virtual ~fractionCalculation() noexcept =0;
};
}  // namespace schemi

#endif /* FRACTIONCALCULATION_HPP_ */
