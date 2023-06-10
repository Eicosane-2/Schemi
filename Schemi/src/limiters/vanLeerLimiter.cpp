/*
 * vanLeerLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vanLeerLimiter.hpp"

schemi::scalar schemi::vanLeerLimiter::vanLeerLimiterCalculation(const scalar r,
		const scalar xiR) const noexcept
{
	scalar vanLeer;

	if (r < 0.)
		vanLeer = 0;
	else
		vanLeer = std::min(2 * r / (1 + r), xiR);

	return vanLeer;
}

schemi::scalar schemi::vanLeerLimiter::vanLeerLimiterCalculation(
		const scalar r) const noexcept
{
	scalar vanLeer;

	if (r < 0.)
		vanLeer = 0;
	else
		vanLeer = 2 * r / (1 + r);

	return vanLeer;
}

schemi::vector schemi::vanLeerLimiter::calculate(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanLeer, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]));

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanLeer.v_r()[j] = vanLeerLimiterCalculation(r.v()[j], xiR.v()[j]);

	return vector { vanLeer.v()[0] * gradientC.v()[0], vanLeer.v()[1]
			* gradientC.v()[1], vanLeer.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::vanLeerLimiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanLeer, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]), 2 / (1 + r.v()[3]), 2 / (1 + r.v()[4]),
			2 / (1 + r.v()[5]), 2 / (1 + r.v()[6]), 2 / (1 + r.v()[7]),
			2 / (1 + r.v()[8]));

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanLeer.v_r()[j] = vanLeerLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor { vanLeer.v()[0] * gradientC.v()[0], vanLeer.v()[1]
			* gradientC.v()[1], vanLeer.v()[2] * gradientC.v()[2],
			vanLeer.v()[3] * gradientC.v()[3], vanLeer.v()[4]
					* gradientC.v()[4], vanLeer.v()[5] * gradientC.v()[5],
			vanLeer.v()[6] * gradientC.v()[6], vanLeer.v()[7]
					* gradientC.v()[7], vanLeer.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::vanLeerLimiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 vanLeer, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
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
		vanLeer.v_r()[j] = vanLeerLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor3 { vanLeer.v()[0] * gradientC.v()[0], vanLeer.v()[1]
			* gradientC.v()[1], vanLeer.v()[2] * gradientC.v()[2],
			vanLeer.v()[3] * gradientC.v()[3], vanLeer.v()[4]
					* gradientC.v()[4], vanLeer.v()[5] * gradientC.v()[5],
			vanLeer.v()[6] * gradientC.v()[6], vanLeer.v()[7]
					* gradientC.v()[7], vanLeer.v()[8] * gradientC.v()[8],
			vanLeer.v()[9] * gradientC.v()[9], vanLeer.v()[10]
					* gradientC.v()[10], vanLeer.v()[11] * gradientC.v()[11],
			vanLeer.v()[12] * gradientC.v()[12], vanLeer.v()[13]
					* gradientC.v()[13], vanLeer.v()[14] * gradientC.v()[14],
			vanLeer.v()[15] * gradientC.v()[15], vanLeer.v()[16]
					* gradientC.v()[16], vanLeer.v()[17] * gradientC.v()[17],
			vanLeer.v()[18] * gradientC.v()[18], vanLeer.v()[19]
					* gradientC.v()[19], vanLeer.v()[20] * gradientC.v()[20],
			vanLeer.v()[21] * gradientC.v()[21], vanLeer.v()[22]
					* gradientC.v()[22], vanLeer.v()[23] * gradientC.v()[23],
			vanLeer.v()[24] * gradientC.v()[24], vanLeer.v()[25]
					* gradientC.v()[25], vanLeer.v()[26] * gradientC.v()[26] };
}

schemi::vector schemi::vanLeerLimiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanLeer;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanLeer.v_r()[j] = vanLeerLimiterCalculation(r.v()[j]);

	return vector { vanLeer.v()[0] * gradientC.v()[0], vanLeer.v()[1]
			* gradientC.v()[1], vanLeer.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::vanLeerLimiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanLeer;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanLeer.v_r()[j] = vanLeerLimiterCalculation(r.v()[j]);

	return tensor { vanLeer.v()[0] * gradientC.v()[0], vanLeer.v()[1]
			* gradientC.v()[1], vanLeer.v()[2] * gradientC.v()[2],
			vanLeer.v()[3] * gradientC.v()[3], vanLeer.v()[4]
					* gradientC.v()[4], vanLeer.v()[5] * gradientC.v()[5],
			vanLeer.v()[6] * gradientC.v()[6], vanLeer.v()[7]
					* gradientC.v()[7], vanLeer.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::vanLeerLimiter::calculateNoRightLimit(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 vanLeer;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		vanLeer.v_r()[j] = vanLeerLimiterCalculation(r.v()[j]);

	return tensor3 { vanLeer.v()[0] * gradientC.v()[0], vanLeer.v()[1]
			* gradientC.v()[1], vanLeer.v()[2] * gradientC.v()[2],
			vanLeer.v()[3] * gradientC.v()[3], vanLeer.v()[4]
					* gradientC.v()[4], vanLeer.v()[5] * gradientC.v()[5],
			vanLeer.v()[6] * gradientC.v()[6], vanLeer.v()[7]
					* gradientC.v()[7], vanLeer.v()[8] * gradientC.v()[8],
			vanLeer.v()[9] * gradientC.v()[9], vanLeer.v()[10]
					* gradientC.v()[10], vanLeer.v()[11] * gradientC.v()[11],
			vanLeer.v()[12] * gradientC.v()[12], vanLeer.v()[13]
					* gradientC.v()[13], vanLeer.v()[14] * gradientC.v()[14],
			vanLeer.v()[15] * gradientC.v()[15], vanLeer.v()[16]
					* gradientC.v()[16], vanLeer.v()[17] * gradientC.v()[17],
			vanLeer.v()[18] * gradientC.v()[18], vanLeer.v()[19]
					* gradientC.v()[19], vanLeer.v()[20] * gradientC.v()[20],
			vanLeer.v()[21] * gradientC.v()[21], vanLeer.v()[22]
					* gradientC.v()[22], vanLeer.v()[23] * gradientC.v()[23],
			vanLeer.v()[24] * gradientC.v()[24], vanLeer.v()[25]
					* gradientC.v()[25], vanLeer.v()[26] * gradientC.v()[26] };
}
