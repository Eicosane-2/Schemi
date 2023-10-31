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
	vector superbee, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	for (std::size_t j = 0; j < vector::vsize; ++j)
		superbee.r()[j] = superbeeLimiterCalculation(r()[j], xiR()[j]);

	return vector { std::get<0>(superbee()) * std::get<0>(gradientC()),
			std::get<1>(superbee()) * std::get<1>(gradientC()), std::get<2>(
					superbee()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::superbeeLimiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor superbee, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		superbee.r()[j] = superbeeLimiterCalculation(r()[j], xiR()[j]);

	return tensor { std::get<0>(superbee()) * std::get<0>(gradientC()),
			std::get<1>(superbee()) * std::get<1>(gradientC()), std::get<2>(
					superbee()) * std::get<2>(gradientC()), std::get<3>(
					superbee()) * std::get<3>(gradientC()), std::get<4>(
					superbee()) * std::get<4>(gradientC()), std::get<5>(
					superbee()) * std::get<5>(gradientC()), std::get<6>(
					superbee()) * std::get<6>(gradientC()), std::get<7>(
					superbee()) * std::get<7>(gradientC()), std::get<8>(
					superbee()) * std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::superbeeLimiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 superbee, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())), 2
			/ (1 + std::get<9>(r())), 2 / (1 + std::get<10>(r())), 2
			/ (1 + std::get<11>(r())), 2 / (1 + std::get<12>(r())), 2
			/ (1 + std::get<13>(r())), 2 / (1 + std::get<14>(r())), 2
			/ (1 + std::get<15>(r())), 2 / (1 + std::get<16>(r())), 2
			/ (1 + std::get<17>(r())), 2 / (1 + std::get<18>(r())), 2
			/ (1 + std::get<19>(r())), 2 / (1 + std::get<20>(r())), 2
			/ (1 + std::get<21>(r())), 2 / (1 + std::get<22>(r())), 2
			/ (1 + std::get<23>(r())), 2 / (1 + std::get<24>(r())), 2
			/ (1 + std::get<25>(r())), 2 / (1 + std::get<26>(r())) };

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		superbee.r()[j] = superbeeLimiterCalculation(r()[j], xiR()[j]);

	return tensor3 { std::get<0>(superbee()) * std::get<0>(gradientC()),
			std::get<1>(superbee()) * std::get<1>(gradientC()), std::get<2>(
					superbee()) * std::get<2>(gradientC()), std::get<3>(
					superbee()) * std::get<3>(gradientC()), std::get<4>(
					superbee()) * std::get<4>(gradientC()), std::get<5>(
					superbee()) * std::get<5>(gradientC()), std::get<6>(
					superbee()) * std::get<6>(gradientC()), std::get<7>(
					superbee()) * std::get<7>(gradientC()), std::get<8>(
					superbee()) * std::get<8>(gradientC()), std::get<9>(
					superbee()) * std::get<9>(gradientC()), std::get<10>(
					superbee()) * std::get<10>(gradientC()), std::get<11>(
					superbee()) * std::get<11>(gradientC()), std::get<12>(
					superbee()) * std::get<12>(gradientC()), std::get<13>(
					superbee()) * std::get<13>(gradientC()), std::get<14>(
					superbee()) * std::get<14>(gradientC()), std::get<15>(
					superbee()) * std::get<15>(gradientC()), std::get<16>(
					superbee()) * std::get<16>(gradientC()), std::get<17>(
					superbee()) * std::get<17>(gradientC()), std::get<18>(
					superbee()) * std::get<18>(gradientC()), std::get<19>(
					superbee()) * std::get<19>(gradientC()), std::get<20>(
					superbee()) * std::get<20>(gradientC()), std::get<21>(
					superbee()) * std::get<21>(gradientC()), std::get<22>(
					superbee()) * std::get<22>(gradientC()), std::get<23>(
					superbee()) * std::get<23>(gradientC()), std::get<24>(
					superbee()) * std::get<24>(gradientC()), std::get<25>(
					superbee()) * std::get<25>(gradientC()), std::get<26>(
					superbee()) * std::get<26>(gradientC()) };
}

schemi::vector schemi::superbeeLimiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector superbee;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		superbee.r()[j] = superbeeLimiterCalculation(r()[j]);

	return vector { std::get<0>(superbee()) * std::get<0>(gradientC()),
			std::get<1>(superbee()) * std::get<1>(gradientC()), std::get<2>(
					superbee()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::superbeeLimiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor superbee;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		superbee.r()[j] = superbeeLimiterCalculation(r()[j]);

	return tensor { std::get<0>(superbee()) * std::get<0>(gradientC()),
			std::get<1>(superbee()) * std::get<1>(gradientC()), std::get<2>(
					superbee()) * std::get<2>(gradientC()), std::get<3>(
					superbee()) * std::get<3>(gradientC()), std::get<4>(
					superbee()) * std::get<4>(gradientC()), std::get<5>(
					superbee()) * std::get<5>(gradientC()), std::get<6>(
					superbee()) * std::get<6>(gradientC()), std::get<7>(
					superbee()) * std::get<7>(gradientC()), std::get<8>(
					superbee()) * std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::superbeeLimiter::calculateNoRightLimit(
		const tensor3 & r, const tensor3 & gradientC) const noexcept
{
	tensor3 superbee;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		superbee.r()[j] = superbeeLimiterCalculation(r()[j]);

	return tensor3 { std::get<0>(superbee()) * std::get<0>(gradientC()),
			std::get<1>(superbee()) * std::get<1>(gradientC()), std::get<2>(
					superbee()) * std::get<2>(gradientC()), std::get<3>(
					superbee()) * std::get<3>(gradientC()), std::get<4>(
					superbee()) * std::get<4>(gradientC()), std::get<5>(
					superbee()) * std::get<5>(gradientC()), std::get<6>(
					superbee()) * std::get<6>(gradientC()), std::get<7>(
					superbee()) * std::get<7>(gradientC()), std::get<8>(
					superbee()) * std::get<8>(gradientC()), std::get<9>(
					superbee()) * std::get<9>(gradientC()), std::get<10>(
					superbee()) * std::get<10>(gradientC()), std::get<11>(
					superbee()) * std::get<11>(gradientC()), std::get<12>(
					superbee()) * std::get<12>(gradientC()), std::get<13>(
					superbee()) * std::get<13>(gradientC()), std::get<14>(
					superbee()) * std::get<14>(gradientC()), std::get<15>(
					superbee()) * std::get<15>(gradientC()), std::get<16>(
					superbee()) * std::get<16>(gradientC()), std::get<17>(
					superbee()) * std::get<17>(gradientC()), std::get<18>(
					superbee()) * std::get<18>(gradientC()), std::get<19>(
					superbee()) * std::get<19>(gradientC()), std::get<20>(
					superbee()) * std::get<20>(gradientC()), std::get<21>(
					superbee()) * std::get<21>(gradientC()), std::get<22>(
					superbee()) * std::get<22>(gradientC()), std::get<23>(
					superbee()) * std::get<23>(gradientC()), std::get<24>(
					superbee()) * std::get<24>(gradientC()), std::get<25>(
					superbee()) * std::get<25>(gradientC()), std::get<26>(
					superbee()) * std::get<26>(gradientC()) };
}
