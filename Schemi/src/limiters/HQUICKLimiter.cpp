/*
 * HQUICKLimiter.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HQUICKLimiter.hpp"

schemi::scalar schemi::HQUICKLimiter::HQUICKLimiterCalculation(const scalar r,
		const scalar xiR) const noexcept
{
	scalar HQUICK;

	if (r <= 0.)
		HQUICK = 0;
	else
		HQUICK = std::min(4 * r / (3 + r), xiR);

	return HQUICK;
}

schemi::scalar schemi::HQUICKLimiter::HQUICKLimiterCalculation(
		const scalar r) const noexcept
{
	scalar HQUICK;

	if (r <= 0.)
		HQUICK = 0;
	else
		HQUICK = 4 * r / (3 + r);

	return HQUICK;
}

schemi::vector schemi::HQUICKLimiter::calculate(const vector & r,
		const vector & gradientC) const noexcept
{
	vector HQUICK, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]));

	for (std::size_t j = 0; j < vector::vsize; ++j)
		HQUICK.v_r()[j] = HQUICKLimiterCalculation(r.v()[j], xiR.v()[j]);

	return vector { HQUICK.v()[0] * gradientC.v()[0], HQUICK.v()[1]
			* gradientC.v()[1], HQUICK.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::HQUICKLimiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor HQUICK, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]), 2 / (1 + r.v()[3]), 2 / (1 + r.v()[4]),
			2 / (1 + r.v()[5]), 2 / (1 + r.v()[6]), 2 / (1 + r.v()[7]),
			2 / (1 + r.v()[8]));

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		HQUICK.v_r()[j] = HQUICKLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor { HQUICK.v()[0] * gradientC.v()[0], HQUICK.v()[1]
			* gradientC.v()[1], HQUICK.v()[2] * gradientC.v()[2], HQUICK.v()[3]
			* gradientC.v()[3], HQUICK.v()[4] * gradientC.v()[4], HQUICK.v()[5]
			* gradientC.v()[5], HQUICK.v()[6] * gradientC.v()[6], HQUICK.v()[7]
			* gradientC.v()[7], HQUICK.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::HQUICKLimiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 HQUICK, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
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
		HQUICK.v_r()[j] = HQUICKLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor3 { HQUICK.v()[0] * gradientC.v()[0], HQUICK.v()[1]
			* gradientC.v()[1], HQUICK.v()[2] * gradientC.v()[2], HQUICK.v()[3]
			* gradientC.v()[3], HQUICK.v()[4] * gradientC.v()[4], HQUICK.v()[5]
			* gradientC.v()[5], HQUICK.v()[6] * gradientC.v()[6], HQUICK.v()[7]
			* gradientC.v()[7], HQUICK.v()[8] * gradientC.v()[8], HQUICK.v()[9]
			* gradientC.v()[9], HQUICK.v()[10] * gradientC.v()[10],
			HQUICK.v()[11] * gradientC.v()[11], HQUICK.v()[12]
					* gradientC.v()[12], HQUICK.v()[13] * gradientC.v()[13],
			HQUICK.v()[14] * gradientC.v()[14], HQUICK.v()[15]
					* gradientC.v()[15], HQUICK.v()[16] * gradientC.v()[16],
			HQUICK.v()[17] * gradientC.v()[17], HQUICK.v()[18]
					* gradientC.v()[18], HQUICK.v()[19] * gradientC.v()[19],
			HQUICK.v()[20] * gradientC.v()[20], HQUICK.v()[21]
					* gradientC.v()[21], HQUICK.v()[22] * gradientC.v()[22],
			HQUICK.v()[23] * gradientC.v()[23], HQUICK.v()[24]
					* gradientC.v()[24], HQUICK.v()[25] * gradientC.v()[25],
			HQUICK.v()[26] * gradientC.v()[26] };
}

schemi::vector schemi::HQUICKLimiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector HQUICK;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		HQUICK.v_r()[j] = HQUICKLimiterCalculation(r.v()[j]);

	return vector { HQUICK.v()[0] * gradientC.v()[0], HQUICK.v()[1]
			* gradientC.v()[1], HQUICK.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::HQUICKLimiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor HQUICK;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		HQUICK.v_r()[j] = HQUICKLimiterCalculation(r.v()[j]);

	return tensor { HQUICK.v()[0] * gradientC.v()[0], HQUICK.v()[1]
			* gradientC.v()[1], HQUICK.v()[2] * gradientC.v()[2], HQUICK.v()[3]
			* gradientC.v()[3], HQUICK.v()[4] * gradientC.v()[4], HQUICK.v()[5]
			* gradientC.v()[5], HQUICK.v()[6] * gradientC.v()[6], HQUICK.v()[7]
			* gradientC.v()[7], HQUICK.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::HQUICKLimiter::calculateNoRightLimit(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 HQUICK;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		HQUICK.v_r()[j] = HQUICKLimiterCalculation(r.v()[j]);

	return tensor3 { HQUICK.v()[0] * gradientC.v()[0], HQUICK.v()[1]
			* gradientC.v()[1], HQUICK.v()[2] * gradientC.v()[2], HQUICK.v()[3]
			* gradientC.v()[3], HQUICK.v()[4] * gradientC.v()[4], HQUICK.v()[5]
			* gradientC.v()[5], HQUICK.v()[6] * gradientC.v()[6], HQUICK.v()[7]
			* gradientC.v()[7], HQUICK.v()[8] * gradientC.v()[8], HQUICK.v()[9]
			* gradientC.v()[9], HQUICK.v()[10] * gradientC.v()[10],
			HQUICK.v()[11] * gradientC.v()[11], HQUICK.v()[12]
					* gradientC.v()[12], HQUICK.v()[13] * gradientC.v()[13],
			HQUICK.v()[14] * gradientC.v()[14], HQUICK.v()[15]
					* gradientC.v()[15], HQUICK.v()[16] * gradientC.v()[16],
			HQUICK.v()[17] * gradientC.v()[17], HQUICK.v()[18]
					* gradientC.v()[18], HQUICK.v()[19] * gradientC.v()[19],
			HQUICK.v()[20] * gradientC.v()[20], HQUICK.v()[21]
					* gradientC.v()[21], HQUICK.v()[22] * gradientC.v()[22],
			HQUICK.v()[23] * gradientC.v()[23], HQUICK.v()[24]
					* gradientC.v()[24], HQUICK.v()[25] * gradientC.v()[25],
			HQUICK.v()[26] * gradientC.v()[26] };
}
