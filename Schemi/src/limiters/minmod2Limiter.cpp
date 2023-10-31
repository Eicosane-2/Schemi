/*
 * minmod2Limiter.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "minmod2Limiter.hpp"

schemi::scalar schemi::minmod2Limiter::minmod2LimiterCalculation(const scalar r,
		const scalar xiR) const noexcept
{
	scalar minmod2;

	if (r <= 0.)
		minmod2 = 0;
	else
		minmod2 = std::min(static_cast<scalar>(1.),
				std::min(4 * r / (1 + r), 2 * xiR));

	return minmod2;
}

schemi::scalar schemi::minmod2Limiter::minmod2LimiterCalculation(
		const scalar r) const noexcept
{
	scalar minmod2;

	if (r <= 0.)
		minmod2 = 0;
	else
		minmod2 = std::min(static_cast<scalar>(1.), (4 * r / (1 + r)));

	return minmod2;
}

schemi::vector schemi::minmod2Limiter::calculate(const vector & r,
		const vector & gradientC) const noexcept
{
	vector minmod2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	for (std::size_t j = 0; j < vector::vsize; ++j)
		minmod2.r()[j] = minmod2LimiterCalculation(r()[j], xiR()[j]);

	return vector { std::get<0>(minmod2()) * std::get<0>(gradientC()), std::get<
			1>(minmod2()) * std::get<1>(gradientC()), std::get<2>(minmod2())
			* std::get<2>(gradientC()) };
}

schemi::tensor schemi::minmod2Limiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor minmod2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		minmod2.r()[j] = minmod2LimiterCalculation(r()[j], xiR()[j]);

	return tensor { std::get<0>(minmod2()) * std::get<0>(gradientC()), std::get<
			1>(minmod2()) * std::get<1>(gradientC()), std::get<2>(minmod2())
			* std::get<2>(gradientC()), std::get<3>(minmod2())
			* std::get<3>(gradientC()), std::get<4>(minmod2())
			* std::get<4>(gradientC()), std::get<5>(minmod2())
			* std::get<5>(gradientC()), std::get<6>(minmod2())
			* std::get<6>(gradientC()), std::get<7>(minmod2())
			* std::get<7>(gradientC()), std::get<8>(minmod2())
			* std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::minmod2Limiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 minmod2, xiR { 2 / (1 + std::get<0>(r())), 2
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
		minmod2.r()[j] = minmod2LimiterCalculation(r()[j], xiR()[j]);

	return tensor3 { std::get<0>(minmod2()) * std::get<0>(gradientC()),
			std::get<1>(minmod2()) * std::get<1>(gradientC()), std::get<2>(
					minmod2()) * std::get<2>(gradientC()), std::get<3>(
					minmod2()) * std::get<3>(gradientC()), std::get<4>(
					minmod2()) * std::get<4>(gradientC()), std::get<5>(
					minmod2()) * std::get<5>(gradientC()), std::get<6>(
					minmod2()) * std::get<6>(gradientC()), std::get<7>(
					minmod2()) * std::get<7>(gradientC()), std::get<8>(
					minmod2()) * std::get<8>(gradientC()), std::get<9>(
					minmod2()) * std::get<9>(gradientC()), std::get<10>(
					minmod2()) * std::get<10>(gradientC()), std::get<11>(
					minmod2()) * std::get<11>(gradientC()), std::get<12>(
					minmod2()) * std::get<12>(gradientC()), std::get<13>(
					minmod2()) * std::get<13>(gradientC()), std::get<14>(
					minmod2()) * std::get<14>(gradientC()), std::get<15>(
					minmod2()) * std::get<15>(gradientC()), std::get<16>(
					minmod2()) * std::get<16>(gradientC()), std::get<17>(
					minmod2()) * std::get<17>(gradientC()), std::get<18>(
					minmod2()) * std::get<18>(gradientC()), std::get<19>(
					minmod2()) * std::get<19>(gradientC()), std::get<20>(
					minmod2()) * std::get<20>(gradientC()), std::get<21>(
					minmod2()) * std::get<21>(gradientC()), std::get<22>(
					minmod2()) * std::get<22>(gradientC()), std::get<23>(
					minmod2()) * std::get<23>(gradientC()), std::get<24>(
					minmod2()) * std::get<24>(gradientC()), std::get<25>(
					minmod2()) * std::get<25>(gradientC()), std::get<26>(
					minmod2()) * std::get<26>(gradientC()) };
}

schemi::vector schemi::minmod2Limiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector minmod2;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		minmod2.r()[j] = minmod2LimiterCalculation(r()[j]);

	return vector { std::get<0>(minmod2()) * std::get<0>(gradientC()), std::get<
			1>(minmod2()) * std::get<1>(gradientC()), std::get<2>(minmod2())
			* std::get<2>(gradientC()) };
}

schemi::tensor schemi::minmod2Limiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor minmod2;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		minmod2.r()[j] = minmod2LimiterCalculation(r()[j]);

	return tensor { std::get<0>(minmod2()) * std::get<0>(gradientC()), std::get<
			1>(minmod2()) * std::get<1>(gradientC()), std::get<2>(minmod2())
			* std::get<2>(gradientC()), std::get<3>(minmod2())
			* std::get<3>(gradientC()), std::get<4>(minmod2())
			* std::get<4>(gradientC()), std::get<5>(minmod2())
			* std::get<5>(gradientC()), std::get<6>(minmod2())
			* std::get<6>(gradientC()), std::get<7>(minmod2())
			* std::get<7>(gradientC()), std::get<8>(minmod2())
			* std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::minmod2Limiter::calculateNoRightLimit(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 minmod2;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		minmod2.r()[j] = minmod2LimiterCalculation(r()[j]);

	return tensor3 { std::get<0>(minmod2()) * std::get<0>(gradientC()),
			std::get<1>(minmod2()) * std::get<1>(gradientC()), std::get<2>(
					minmod2()) * std::get<2>(gradientC()), std::get<3>(
					minmod2()) * std::get<3>(gradientC()), std::get<4>(
					minmod2()) * std::get<4>(gradientC()), std::get<5>(
					minmod2()) * std::get<5>(gradientC()), std::get<6>(
					minmod2()) * std::get<6>(gradientC()), std::get<7>(
					minmod2()) * std::get<7>(gradientC()), std::get<8>(
					minmod2()) * std::get<8>(gradientC()), std::get<9>(
					minmod2()) * std::get<9>(gradientC()), std::get<10>(
					minmod2()) * std::get<10>(gradientC()), std::get<11>(
					minmod2()) * std::get<11>(gradientC()), std::get<12>(
					minmod2()) * std::get<12>(gradientC()), std::get<13>(
					minmod2()) * std::get<13>(gradientC()), std::get<14>(
					minmod2()) * std::get<14>(gradientC()), std::get<15>(
					minmod2()) * std::get<15>(gradientC()), std::get<16>(
					minmod2()) * std::get<16>(gradientC()), std::get<17>(
					minmod2()) * std::get<17>(gradientC()), std::get<18>(
					minmod2()) * std::get<18>(gradientC()), std::get<19>(
					minmod2()) * std::get<19>(gradientC()), std::get<20>(
					minmod2()) * std::get<20>(gradientC()), std::get<21>(
					minmod2()) * std::get<21>(gradientC()), std::get<22>(
					minmod2()) * std::get<22>(gradientC()), std::get<23>(
					minmod2()) * std::get<23>(gradientC()), std::get<24>(
					minmod2()) * std::get<24>(gradientC()), std::get<25>(
					minmod2()) * std::get<25>(gradientC()), std::get<26>(
					minmod2()) * std::get<26>(gradientC()) };
}
