/*
 * fractionCalculation.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "fractionCalculation.hpp"

std::valarray<std::valarray<schemi::scalar>> schemi::fractionCalculation::calcMolarFrac(
		const std::vector<const std::valarray<scalar>*> & concentrations) const noexcept
{
	std::valarray<std::valarray<scalar>> XOutput { concentrations.size() - 1 };

	for (std::size_t k = 0; k < XOutput.size(); ++k)
		XOutput[k] = (*concentrations[k + 1]) / (*concentrations[0]);

	return XOutput;
}

std::valarray<std::valarray<schemi::scalar>> schemi::fractionCalculation::rearrangeMolFrac(
		const std::valarray<std::valarray<scalar>> & in) const noexcept
{
	std::valarray<std::valarray<scalar>> out { in[0].size() };

	for (std::size_t i = 0; i < out.size(); ++i)
	{
		out[i].resize(in.size());

		for (std::size_t k = 0; k < in.size(); ++k)
			out[i][k] = in[k][i];
	}

	return out;
}

std::valarray<schemi::scalar> schemi::fractionCalculation::calcMolarFrac(
		const std::valarray<scalar> & concentrations) const noexcept
{
	std::valarray<scalar> XOutput(concentrations.size() - 1);

	for (std::size_t k = 0; k < XOutput.size(); ++k)
		XOutput[k] = concentrations[k + 1] / concentrations[0];

	return XOutput;
}

schemi::fractionCalculation::~fractionCalculation() noexcept
{
}
