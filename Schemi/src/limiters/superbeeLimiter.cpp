/*
 * superbeeLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "superbeeLimiter.hpp"

schemi::scalar schemi::superbeeLimiter::superbeeLimiterCalculation(
		const scalar r, const scalar xiR) const noexcept
{
	scalar superbee;

	if (r <= 0.)
		superbee = 0;
	else if (r <= 0.5)
		superbee = 2 * r;
	else if (r <= 1.)
		superbee = r;
	else
		superbee = std::min(std::min(r, static_cast<scalar>(2.)), xiR);

	return superbee;
}

schemi::scalar schemi::superbeeLimiter::superbeeLimiterCalculation(
		const scalar r) const noexcept
{
	scalar superbee;

	if (r <= 0.)
		superbee = 0;
	else if (r <= 0.5)
		superbee = 2 * r;
	else if (r <= 1.)
		superbee = r;
	else
		superbee = std::min(r, static_cast<scalar>(2.));

	return superbee;
}

schemi::vector schemi::superbeeLimiter::calculate(const vector & r,
		const vector & gradientC) const noexcept
{
	vector superbee, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]));

	for (std::size_t j = 0; j < vector::vsize; ++j)
		superbee.v_r()[j] = superbeeLimiterCalculation(r.v()[j], xiR.v()[j]);

	return vector { superbee.v()[0] * gradientC.v()[0], superbee.v()[1]
			* gradientC.v()[1], superbee.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::superbeeLimiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor superbee, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]), 2 / (1 + r.v()[3]), 2 / (1 + r.v()[4]),
			2 / (1 + r.v()[5]), 2 / (1 + r.v()[6]), 2 / (1 + r.v()[7]),
			2 / (1 + r.v()[8]));

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		superbee.v_r()[j] = superbeeLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor { superbee.v()[0] * gradientC.v()[0], superbee.v()[1]
			* gradientC.v()[1], superbee.v()[2] * gradientC.v()[2],
			superbee.v()[3] * gradientC.v()[3], superbee.v()[4]
					* gradientC.v()[4], superbee.v()[5] * gradientC.v()[5],
			superbee.v()[6] * gradientC.v()[6], superbee.v()[7]
					* gradientC.v()[7], superbee.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::superbeeLimiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 superbee, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
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
		superbee.v_r()[j] = superbeeLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor3 { superbee.v()[0] * gradientC.v()[0], superbee.v()[1]
			* gradientC.v()[1], superbee.v()[2] * gradientC.v()[2],
			superbee.v()[3] * gradientC.v()[3], superbee.v()[4]
					* gradientC.v()[4], superbee.v()[5] * gradientC.v()[5],
			superbee.v()[6] * gradientC.v()[6], superbee.v()[7]
					* gradientC.v()[7], superbee.v()[8] * gradientC.v()[8],
			superbee.v()[9] * gradientC.v()[9], superbee.v()[10]
					* gradientC.v()[10], superbee.v()[11] * gradientC.v()[11],
			superbee.v()[12] * gradientC.v()[12], superbee.v()[13]
					* gradientC.v()[13], superbee.v()[14] * gradientC.v()[14],
			superbee.v()[15] * gradientC.v()[15], superbee.v()[16]
					* gradientC.v()[16], superbee.v()[17] * gradientC.v()[17],
			superbee.v()[18] * gradientC.v()[18], superbee.v()[19]
					* gradientC.v()[19], superbee.v()[20] * gradientC.v()[20],
			superbee.v()[21] * gradientC.v()[21], superbee.v()[22]
					* gradientC.v()[22], superbee.v()[23] * gradientC.v()[23],
			superbee.v()[24] * gradientC.v()[24], superbee.v()[25]
					* gradientC.v()[25], superbee.v()[26] * gradientC.v()[26] };
}

schemi::vector schemi::superbeeLimiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector superbee;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		superbee.v_r()[j] = superbeeLimiterCalculation(r.v()[j]);

	return vector { superbee.v()[0] * gradientC.v()[0], superbee.v()[1]
			* gradientC.v()[1], superbee.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::superbeeLimiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor superbee;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		superbee.v_r()[j] = superbeeLimiterCalculation(r.v()[j]);

	return tensor { superbee.v()[0] * gradientC.v()[0], superbee.v()[1]
			* gradientC.v()[1], superbee.v()[2] * gradientC.v()[2],
			superbee.v()[3] * gradientC.v()[3], superbee.v()[4]
					* gradientC.v()[4], superbee.v()[5] * gradientC.v()[5],
			superbee.v()[6] * gradientC.v()[6], superbee.v()[7]
					* gradientC.v()[7], superbee.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::superbeeLimiter::calculateNoRightLimit(
		const tensor3 & r, const tensor3 & gradientC) const noexcept
{
	tensor3 superbee;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		superbee.v_r()[j] = superbeeLimiterCalculation(r.v()[j]);

	return tensor3 { superbee.v()[0] * gradientC.v()[0], superbee.v()[1]
			* gradientC.v()[1], superbee.v()[2] * gradientC.v()[2],
			superbee.v()[3] * gradientC.v()[3], superbee.v()[4]
					* gradientC.v()[4], superbee.v()[5] * gradientC.v()[5],
			superbee.v()[6] * gradientC.v()[6], superbee.v()[7]
					* gradientC.v()[7], superbee.v()[8] * gradientC.v()[8],
			superbee.v()[9] * gradientC.v()[9], superbee.v()[10]
					* gradientC.v()[10], superbee.v()[11] * gradientC.v()[11],
			superbee.v()[12] * gradientC.v()[12], superbee.v()[13]
					* gradientC.v()[13], superbee.v()[14] * gradientC.v()[14],
			superbee.v()[15] * gradientC.v()[15], superbee.v()[16]
					* gradientC.v()[16], superbee.v()[17] * gradientC.v()[17],
			superbee.v()[18] * gradientC.v()[18], superbee.v()[19]
					* gradientC.v()[19], superbee.v()[20] * gradientC.v()[20],
			superbee.v()[21] * gradientC.v()[21], superbee.v()[22]
					* gradientC.v()[22], superbee.v()[23] * gradientC.v()[23],
			superbee.v()[24] * gradientC.v()[24], superbee.v()[25]
					* gradientC.v()[25], superbee.v()[26] * gradientC.v()[26] };
}
