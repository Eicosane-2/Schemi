/*
 * vanLeer2Limiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vanLeer2Limiter.hpp"

#include "intExpPow.hpp"

schemi::scalar schemi::vanLeer2Limiter::vanLeer2LimiterCalculation(
		const scalar r, const scalar xiR) const noexcept
{
	scalar vanLeer2;

	if (r <= 0.)
		vanLeer2 = 0;
	else
	{
		const auto denom { 1 + r };

		vanLeer2 = std::min(4 * r / pow<scalar, 2>(denom), 2 * xiR);
	}

	return vanLeer2;
}

schemi::scalar schemi::vanLeer2Limiter::vanLeer2LimiterCalculation(
		const scalar r) const noexcept
{
	scalar vanLeer2;

	if (r <= 0.)
		vanLeer2 = 0;
	else
	{
		const auto denom { 1 + r };

		vanLeer2 = 4 * r / pow<scalar, 2>(denom);
	}

	return vanLeer2;
}

schemi::vector schemi::vanLeer2Limiter::calculate(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanLeer2, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]));

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanLeer2.v_r()[j] = vanLeer2LimiterCalculation(r.v()[j], xiR.v()[j]);

	return vector { vanLeer2.v()[0] * gradientC.v()[0], vanLeer2.v()[1]
			* gradientC.v()[1], vanLeer2.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::vanLeer2Limiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanLeer2, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]), 2 / (1 + r.v()[3]), 2 / (1 + r.v()[4]),
			2 / (1 + r.v()[5]), 2 / (1 + r.v()[6]), 2 / (1 + r.v()[7]),
			2 / (1 + r.v()[8]));

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanLeer2.v_r()[j] = vanLeer2LimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor { vanLeer2.v()[0] * gradientC.v()[0], vanLeer2.v()[1]
			* gradientC.v()[1], vanLeer2.v()[2] * gradientC.v()[2],
			vanLeer2.v()[3] * gradientC.v()[3], vanLeer2.v()[4]
					* gradientC.v()[4], vanLeer2.v()[5] * gradientC.v()[5],
			vanLeer2.v()[6] * gradientC.v()[6], vanLeer2.v()[7]
					* gradientC.v()[7], vanLeer2.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::vanLeer2Limiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 vanLeer2, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]), 2 / (1 + r.v()[3]), 2 / (1 + r.v()[4]),
			2 / (1 + r.v()[5]), 2 / (1 + r.v()[6]), 2 / (1 + r.v()[7]),
			2 / (1 + r.v()[8]), 2 / (1 + r.v()[9]), 2 / (1 + r.v()[10]),
			2 / (1 + r.v()[11]), 2 / (1 + r.v()[12]), 2 / (1 + r.v()[13]),
			2 / (1 + r.v()[14]), 2 / (1 + r.v()[15]), 2 / (1 + r.v()[16]),
			2 / (1 + r.v()[17]), 2 / (1 + r.v()[18]), 2 / (1 + r.v()[19]),
			2 / (1 + r.v()[20]), 2 / (1 + r.v()[21]), 2 / (1 + r.v()[22]),
			2 / (1 + r.v()[23]), 2 / (1 + r.v()[24]), 2 / (1 + r.v()[25]),
			2 / (1 + r.v()[26]));

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		vanLeer2.v_r()[j] = vanLeer2LimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor3 { vanLeer2.v()[0] * gradientC.v()[0], vanLeer2.v()[1]
			* gradientC.v()[1], vanLeer2.v()[2] * gradientC.v()[2],
			vanLeer2.v()[3] * gradientC.v()[3], vanLeer2.v()[4]
					* gradientC.v()[4], vanLeer2.v()[5] * gradientC.v()[5],
			vanLeer2.v()[6] * gradientC.v()[6], vanLeer2.v()[7]
					* gradientC.v()[7], vanLeer2.v()[8] * gradientC.v()[8],
			vanLeer2.v()[9] * gradientC.v()[9], vanLeer2.v()[10]
					* gradientC.v()[10], vanLeer2.v()[11] * gradientC.v()[11],
			vanLeer2.v()[12] * gradientC.v()[12], vanLeer2.v()[13]
					* gradientC.v()[13], vanLeer2.v()[14] * gradientC.v()[14],
			vanLeer2.v()[15] * gradientC.v()[15], vanLeer2.v()[16]
					* gradientC.v()[16], vanLeer2.v()[17] * gradientC.v()[17],
			vanLeer2.v()[18] * gradientC.v()[18], vanLeer2.v()[19]
					* gradientC.v()[19], vanLeer2.v()[20] * gradientC.v()[20],
			vanLeer2.v()[21] * gradientC.v()[21], vanLeer2.v()[22]
					* gradientC.v()[22], vanLeer2.v()[23] * gradientC.v()[23],
			vanLeer2.v()[24] * gradientC.v()[24], vanLeer2.v()[25]
					* gradientC.v()[25], vanLeer2.v()[26] * gradientC.v()[26] };
}

schemi::vector schemi::vanLeer2Limiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanLeer2;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanLeer2.v_r()[j] = vanLeer2LimiterCalculation(r.v()[j]);

	return vector { vanLeer2.v()[0] * gradientC.v()[0], vanLeer2.v()[1]
			* gradientC.v()[1], vanLeer2.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::vanLeer2Limiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanLeer2;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanLeer2.v_r()[j] = vanLeer2LimiterCalculation(r.v()[j]);

	return tensor { vanLeer2.v()[0] * gradientC.v()[0], vanLeer2.v()[1]
			* gradientC.v()[1], vanLeer2.v()[2] * gradientC.v()[2],
			vanLeer2.v()[3] * gradientC.v()[3], vanLeer2.v()[4]
					* gradientC.v()[4], vanLeer2.v()[5] * gradientC.v()[5],
			vanLeer2.v()[6] * gradientC.v()[6], vanLeer2.v()[7]
					* gradientC.v()[7], vanLeer2.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::vanLeer2Limiter::calculateNoRightLimit(
		const tensor3 & r, const tensor3 & gradientC) const noexcept
{
	tensor3 vanLeer2;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		vanLeer2.v_r()[j] = vanLeer2LimiterCalculation(r.v()[j]);

	return tensor3 { vanLeer2.v()[0] * gradientC.v()[0], vanLeer2.v()[1]
			* gradientC.v()[1], vanLeer2.v()[2] * gradientC.v()[2],
			vanLeer2.v()[3] * gradientC.v()[3], vanLeer2.v()[4]
					* gradientC.v()[4], vanLeer2.v()[5] * gradientC.v()[5],
			vanLeer2.v()[6] * gradientC.v()[6], vanLeer2.v()[7]
					* gradientC.v()[7], vanLeer2.v()[8] * gradientC.v()[8],
			vanLeer2.v()[9] * gradientC.v()[9], vanLeer2.v()[10]
					* gradientC.v()[10], vanLeer2.v()[11] * gradientC.v()[11],
			vanLeer2.v()[12] * gradientC.v()[12], vanLeer2.v()[13]
					* gradientC.v()[13], vanLeer2.v()[14] * gradientC.v()[14],
			vanLeer2.v()[15] * gradientC.v()[15], vanLeer2.v()[16]
					* gradientC.v()[16], vanLeer2.v()[17] * gradientC.v()[17],
			vanLeer2.v()[18] * gradientC.v()[18], vanLeer2.v()[19]
					* gradientC.v()[19], vanLeer2.v()[20] * gradientC.v()[20],
			vanLeer2.v()[21] * gradientC.v()[21], vanLeer2.v()[22]
					* gradientC.v()[22], vanLeer2.v()[23] * gradientC.v()[23],
			vanLeer2.v()[24] * gradientC.v()[24], vanLeer2.v()[25]
					* gradientC.v()[25], vanLeer2.v()[26] * gradientC.v()[26] };
}
