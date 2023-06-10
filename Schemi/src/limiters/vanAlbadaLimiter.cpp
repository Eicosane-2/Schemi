/*
 * vanAlbadaLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vanAlbadaLimiter.hpp"

#include "intExpPow.hpp"

schemi::scalar schemi::vanAlbadaLimiter::vanAlbadaLimiterCalculation(
		const scalar r, const scalar xiR) const noexcept
{
	scalar vanAlbada;

	if (r <= 0.)
		vanAlbada = 0;
	else
		vanAlbada = std::min((pow<scalar, 2>(r) + r) / (1 + pow<scalar, 2>(r)),
				xiR);

	return vanAlbada;
}

schemi::scalar schemi::vanAlbadaLimiter::vanAlbadaLimiterCalculation(
		const scalar r) const noexcept
{
	scalar vanAlbada;

	if (r <= 0.)
		vanAlbada = 0;
	else
		vanAlbada = (pow<scalar, 2>(r) + r) / (1 + pow<scalar, 2>(r));

	return vanAlbada;
}

schemi::vector schemi::vanAlbadaLimiter::calculate(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanAlbada, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]));

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanAlbada.v_r()[j] = vanAlbadaLimiterCalculation(r.v()[j], xiR.v()[j]);

	return vector { vanAlbada.v()[0] * gradientC.v()[0], vanAlbada.v()[1]
			* gradientC.v()[1], vanAlbada.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::vanAlbadaLimiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanAlbada, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
			2 / (1 + r.v()[2]), 2 / (1 + r.v()[3]), 2 / (1 + r.v()[4]),
			2 / (1 + r.v()[5]), 2 / (1 + r.v()[6]), 2 / (1 + r.v()[7]),
			2 / (1 + r.v()[8]));

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanAlbada.v_r()[j] = vanAlbadaLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor { vanAlbada.v()[0] * gradientC.v()[0], vanAlbada.v()[1]
			* gradientC.v()[1], vanAlbada.v()[2] * gradientC.v()[2],
			vanAlbada.v()[3] * gradientC.v()[3], vanAlbada.v()[4]
					* gradientC.v()[4], vanAlbada.v()[5] * gradientC.v()[5],
			vanAlbada.v()[6] * gradientC.v()[6], vanAlbada.v()[7]
					* gradientC.v()[7], vanAlbada.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::vanAlbadaLimiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 vanAlbada, xiR(2 / (1 + r.v()[0]), 2 / (1 + r.v()[1]),
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
		vanAlbada.v_r()[j] = vanAlbadaLimiterCalculation(r.v()[j], xiR.v()[j]);

	return tensor3 { vanAlbada.v()[0] * gradientC.v()[0], vanAlbada.v()[1]
			* gradientC.v()[1], vanAlbada.v()[2] * gradientC.v()[2],
			vanAlbada.v()[3] * gradientC.v()[3], vanAlbada.v()[4]
					* gradientC.v()[4], vanAlbada.v()[5] * gradientC.v()[5],
			vanAlbada.v()[6] * gradientC.v()[6], vanAlbada.v()[7]
					* gradientC.v()[7], vanAlbada.v()[8] * gradientC.v()[8],
			vanAlbada.v()[9] * gradientC.v()[9], vanAlbada.v()[10]
					* gradientC.v()[10], vanAlbada.v()[11] * gradientC.v()[11],
			vanAlbada.v()[12] * gradientC.v()[12], vanAlbada.v()[13]
					* gradientC.v()[13], vanAlbada.v()[14] * gradientC.v()[14],
			vanAlbada.v()[15] * gradientC.v()[15], vanAlbada.v()[16]
					* gradientC.v()[16], vanAlbada.v()[17] * gradientC.v()[17],
			vanAlbada.v()[18] * gradientC.v()[18], vanAlbada.v()[19]
					* gradientC.v()[19], vanAlbada.v()[20] * gradientC.v()[20],
			vanAlbada.v()[21] * gradientC.v()[21], vanAlbada.v()[22]
					* gradientC.v()[22], vanAlbada.v()[23] * gradientC.v()[23],
			vanAlbada.v()[24] * gradientC.v()[24], vanAlbada.v()[25]
					* gradientC.v()[25], vanAlbada.v()[26] * gradientC.v()[26] };
}

schemi::vector schemi::vanAlbadaLimiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanAlbada;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanAlbada.v_r()[j] = vanAlbadaLimiterCalculation(r.v()[j]);

	return vector { vanAlbada.v()[0] * gradientC.v()[0], vanAlbada.v()[1]
			* gradientC.v()[1], vanAlbada.v()[2] * gradientC.v()[2] };
}

schemi::tensor schemi::vanAlbadaLimiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanAlbada;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanAlbada.v_r()[j] = vanAlbadaLimiterCalculation(r.v()[j]);

	return tensor { vanAlbada.v()[0] * gradientC.v()[0], vanAlbada.v()[1]
			* gradientC.v()[1], vanAlbada.v()[2] * gradientC.v()[2],
			vanAlbada.v()[3] * gradientC.v()[3], vanAlbada.v()[4]
					* gradientC.v()[4], vanAlbada.v()[5] * gradientC.v()[5],
			vanAlbada.v()[6] * gradientC.v()[6], vanAlbada.v()[7]
					* gradientC.v()[7], vanAlbada.v()[8] * gradientC.v()[8] };
}

schemi::tensor3 schemi::vanAlbadaLimiter::calculateNoRightLimit(
		const tensor3 & r, const tensor3 & gradientC) const noexcept
{
	tensor3 vanAlbada;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		vanAlbada.v_r()[j] = vanAlbadaLimiterCalculation(r.v()[j]);

	return tensor3 { vanAlbada.v()[0] * gradientC.v()[0], vanAlbada.v()[1]
			* gradientC.v()[1], vanAlbada.v()[2] * gradientC.v()[2],
			vanAlbada.v()[3] * gradientC.v()[3], vanAlbada.v()[4]
					* gradientC.v()[4], vanAlbada.v()[5] * gradientC.v()[5],
			vanAlbada.v()[6] * gradientC.v()[6], vanAlbada.v()[7]
					* gradientC.v()[7], vanAlbada.v()[8] * gradientC.v()[8],
			vanAlbada.v()[9] * gradientC.v()[9], vanAlbada.v()[10]
					* gradientC.v()[10], vanAlbada.v()[11] * gradientC.v()[11],
			vanAlbada.v()[12] * gradientC.v()[12], vanAlbada.v()[13]
					* gradientC.v()[13], vanAlbada.v()[14] * gradientC.v()[14],
			vanAlbada.v()[15] * gradientC.v()[15], vanAlbada.v()[16]
					* gradientC.v()[16], vanAlbada.v()[17] * gradientC.v()[17],
			vanAlbada.v()[18] * gradientC.v()[18], vanAlbada.v()[19]
					* gradientC.v()[19], vanAlbada.v()[20] * gradientC.v()[20],
			vanAlbada.v()[21] * gradientC.v()[21], vanAlbada.v()[22]
					* gradientC.v()[22], vanAlbada.v()[23] * gradientC.v()[23],
			vanAlbada.v()[24] * gradientC.v()[24], vanAlbada.v()[25]
					* gradientC.v()[25], vanAlbada.v()[26] * gradientC.v()[26] };
}
