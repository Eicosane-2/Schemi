/*
 * minmodLimiter.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "minmodLimiter.hpp"

schemi::scalar schemi::minmodLimiter::minmodLimiterCalculation(const scalar r,
		const scalar xiR) const noexcept
{
	scalar minmod;

	if (r <= 0.)
		minmod = 0;
	else if (r <= 1.)
		minmod = r;
	else
		minmod = std::min(1., xiR);

	return minmod;
}

schemi::scalar schemi::minmodLimiter::minmodLimiterCalculation(
		const scalar r) const noexcept
{
	scalar minmod;

	if (r <= 0.)
		minmod = 0;
	else if (r <= 1.)
		minmod = r;
	else
		minmod = 1.;

	return minmod;
}

schemi::vector schemi::minmodLimiter::calculate(const vector & r,
		const vector & gradientC) const noexcept
{
	vector minmod, xiR { 2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]), 2
			/ (1 + r.v()[2]) };

	for (std::size_t j = 0; j < vector::vsize; ++j)
		minmod.v_r()[j] = minmodLimiterCalculation(r.v()[j], xiR.v()[j]);

	return vector { minmod.v()[0] * gradientC.v()[0], minmod.v()[1]
			* gradientC.v()[1], minmod.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::minmodLimiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor minmod, xiR { 2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]), 2
			/ (1 + r.v()[2]), 2 / (1 + r.v()[3]), 2 / (1 + r.v()[4]), 2
			/ (1 + r.v()[5]), 2 / (1 + r.v()[6]), 2 / (1 + r.v()[7]), 2
			/ (1 + r.v()[8]) };

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		minmod.v_r()[j] = minmodLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor { minmod.v()[0] * gradientC.v()[0], minmod.v()[1]
			* gradientC.v()[1], minmod.v()[2] * gradientC.v()[2], minmod.v()[3]
			* gradientC.v()[3], minmod.v()[4] * gradientC.v()[4], minmod.v()[5]
			* gradientC.v()[5], minmod.v()[6] * gradientC.v()[6], minmod.v()[7]
			* gradientC.v()[7], minmod.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::minmodLimiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 minmod, xiR { 2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]), 2
			/ (1 + r.v()[2]), 2 / (1 + r.v()[3]), 2 / (1 + r.v()[4]), 2
			/ (1 + r.v()[5]), 2 / (1 + r.v()[6]), 2 / (1 + r.v()[7]), 2
			/ (1 + r.v()[8]), 2 / (1 + r.v()[9]), 2 / (1 + r.v()[10]), 2
			/ (1 + r.v()[11]), 2 / (1 + r.v()[12]), 2 / (1 + r.v()[13]), 2
			/ (1 + r.v()[14]), 2 / (1 + r.v()[15]), 2 / (1 + r.v()[16]), 2
			/ (1 + r.v()[17]), 2 / (1 + r.v()[18]), 2 / (1 + r.v()[19]), 2
			/ (1 + r.v()[20]), 2 / (1 + r.v()[21]), 2 / (1 + r.v()[22]), 2
			/ (1 + r.v()[23]), 2 / (1 + r.v()[24]), 2 / (1 + r.v()[25]), 2
			/ (1 + r.v()[26]) };

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		minmod.v_r()[j] = minmodLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor3 { minmod.v()[0] * gradientC.v()[0], minmod.v()[1]
			* gradientC.v()[1], minmod.v()[2] * gradientC.v()[2], minmod.v()[3]
			* gradientC.v()[3], minmod.v()[4] * gradientC.v()[4], minmod.v()[5]
			* gradientC.v()[5], minmod.v()[6] * gradientC.v()[6], minmod.v()[7]
			* gradientC.v()[7], minmod.v()[8] * gradientC.v()[8], minmod.v()[9]
			* gradientC.v()[9], minmod.v()[10] * gradientC.v()[10],
			minmod.v()[11] * gradientC.v()[11], minmod.v()[12]
					* gradientC.v()[12], minmod.v()[13] * gradientC.v()[13],
			minmod.v()[14] * gradientC.v()[14], minmod.v()[15]
					* gradientC.v()[15], minmod.v()[16] * gradientC.v()[16],
			minmod.v()[17] * gradientC.v()[17], minmod.v()[18]
					* gradientC.v()[18], minmod.v()[19] * gradientC.v()[19],
			minmod.v()[20] * gradientC.v()[20], minmod.v()[21]
					* gradientC.v()[21], minmod.v()[22] * gradientC.v()[22],
			minmod.v()[23] * gradientC.v()[23], minmod.v()[24]
					* gradientC.v()[24], minmod.v()[25] * gradientC.v()[25],
			minmod.v()[26] * gradientC.v()[26] };
}

schemi::vector schemi::minmodLimiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector minmod;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		minmod.v_r()[j] = minmodLimiterCalculation(r.v()[j]);

	return vector { minmod.v()[0] * gradientC.v()[0], minmod.v()[1]
			* gradientC.v()[1], minmod.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::minmodLimiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor minmod;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		minmod.v_r()[j] = minmodLimiterCalculation(r.v()[j]);

	return tensor { minmod.v()[0] * gradientC.v()[0], minmod.v()[1]
			* gradientC.v()[1], minmod.v()[2] * gradientC.v()[2], minmod.v()[3]
			* gradientC.v()[3], minmod.v()[4] * gradientC.v()[4], minmod.v()[5]
			* gradientC.v()[5], minmod.v()[6] * gradientC.v()[6], minmod.v()[7]
			* gradientC.v()[7], minmod.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::minmodLimiter::calculateNoRightLimit(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 minmod;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		minmod.v_r()[j] = minmodLimiterCalculation(r.v()[j]);

	return tensor3 { minmod.v()[0] * gradientC.v()[0], minmod.v()[1]
			* gradientC.v()[1], minmod.v()[2] * gradientC.v()[2], minmod.v()[3]
			* gradientC.v()[3], minmod.v()[4] * gradientC.v()[4], minmod.v()[5]
			* gradientC.v()[5], minmod.v()[6] * gradientC.v()[6], minmod.v()[7]
			* gradientC.v()[7], minmod.v()[8] * gradientC.v()[8], minmod.v()[9]
			* gradientC.v()[9], minmod.v()[10] * gradientC.v()[10],
			minmod.v()[11] * gradientC.v()[11], minmod.v()[12]
					* gradientC.v()[12], minmod.v()[13] * gradientC.v()[13],
			minmod.v()[14] * gradientC.v()[14], minmod.v()[15]
					* gradientC.v()[15], minmod.v()[16] * gradientC.v()[16],
			minmod.v()[17] * gradientC.v()[17], minmod.v()[18]
					* gradientC.v()[18], minmod.v()[19] * gradientC.v()[19],
			minmod.v()[20] * gradientC.v()[20], minmod.v()[21]
					* gradientC.v()[21], minmod.v()[22] * gradientC.v()[22],
			minmod.v()[23] * gradientC.v()[23], minmod.v()[24]
					* gradientC.v()[24], minmod.v()[25] * gradientC.v()[25],
			minmod.v()[26] * gradientC.v()[26] };
}
