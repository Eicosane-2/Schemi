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
	vector vanLeer2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanLeer2.r()[j] = vanLeer2LimiterCalculation(r()[j], xiR()[j]);

	return vector { std::get<0>(vanLeer2()) * std::get<0>(gradientC()),
			std::get<1>(vanLeer2()) * std::get<1>(gradientC()), std::get<2>(
					vanLeer2()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::vanLeer2Limiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanLeer2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanLeer2.r()[j] = vanLeer2LimiterCalculation(r()[j], xiR()[j]);

	return tensor { std::get<0>(vanLeer2()) * std::get<0>(gradientC()),
			std::get<1>(vanLeer2()) * std::get<1>(gradientC()), std::get<2>(
					vanLeer2()) * std::get<2>(gradientC()), std::get<3>(
					vanLeer2()) * std::get<3>(gradientC()), std::get<4>(
					vanLeer2()) * std::get<4>(gradientC()), std::get<5>(
					vanLeer2()) * std::get<5>(gradientC()), std::get<6>(
					vanLeer2()) * std::get<6>(gradientC()), std::get<7>(
					vanLeer2()) * std::get<7>(gradientC()), std::get<8>(
					vanLeer2()) * std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::vanLeer2Limiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 vanLeer2, xiR { 2 / (1 + std::get<0>(r())), 2
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
		vanLeer2.r()[j] = vanLeer2LimiterCalculation(r()[j], xiR()[j]);

	return tensor3 { std::get<0>(vanLeer2()) * std::get<0>(gradientC()),
			std::get<1>(vanLeer2()) * std::get<1>(gradientC()), std::get<2>(
					vanLeer2()) * std::get<2>(gradientC()), std::get<3>(
					vanLeer2()) * std::get<3>(gradientC()), std::get<4>(
					vanLeer2()) * std::get<4>(gradientC()), std::get<5>(
					vanLeer2()) * std::get<5>(gradientC()), std::get<6>(
					vanLeer2()) * std::get<6>(gradientC()), std::get<7>(
					vanLeer2()) * std::get<7>(gradientC()), std::get<8>(
					vanLeer2()) * std::get<8>(gradientC()), std::get<9>(
					vanLeer2()) * std::get<9>(gradientC()), std::get<10>(
					vanLeer2()) * std::get<10>(gradientC()), std::get<11>(
					vanLeer2()) * std::get<11>(gradientC()), std::get<12>(
					vanLeer2()) * std::get<12>(gradientC()), std::get<13>(
					vanLeer2()) * std::get<13>(gradientC()), std::get<14>(
					vanLeer2()) * std::get<14>(gradientC()), std::get<15>(
					vanLeer2()) * std::get<15>(gradientC()), std::get<16>(
					vanLeer2()) * std::get<16>(gradientC()), std::get<17>(
					vanLeer2()) * std::get<17>(gradientC()), std::get<18>(
					vanLeer2()) * std::get<18>(gradientC()), std::get<19>(
					vanLeer2()) * std::get<19>(gradientC()), std::get<20>(
					vanLeer2()) * std::get<20>(gradientC()), std::get<21>(
					vanLeer2()) * std::get<21>(gradientC()), std::get<22>(
					vanLeer2()) * std::get<22>(gradientC()), std::get<23>(
					vanLeer2()) * std::get<23>(gradientC()), std::get<24>(
					vanLeer2()) * std::get<24>(gradientC()), std::get<25>(
					vanLeer2()) * std::get<25>(gradientC()), std::get<26>(
					vanLeer2()) * std::get<26>(gradientC()) };
}

schemi::vector schemi::vanLeer2Limiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanLeer2;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanLeer2.r()[j] = vanLeer2LimiterCalculation(r()[j]);

	return vector { std::get<0>(vanLeer2()) * std::get<0>(gradientC()),
			std::get<1>(vanLeer2()) * std::get<1>(gradientC()), std::get<2>(
					vanLeer2()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::vanLeer2Limiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanLeer2;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanLeer2.r()[j] = vanLeer2LimiterCalculation(r()[j]);

	return tensor { std::get<0>(vanLeer2()) * std::get<0>(gradientC()),
			std::get<1>(vanLeer2()) * std::get<1>(gradientC()), std::get<2>(
					vanLeer2()) * std::get<2>(gradientC()), std::get<3>(
					vanLeer2()) * std::get<3>(gradientC()), std::get<4>(
					vanLeer2()) * std::get<4>(gradientC()), std::get<5>(
					vanLeer2()) * std::get<5>(gradientC()), std::get<6>(
					vanLeer2()) * std::get<6>(gradientC()), std::get<7>(
					vanLeer2()) * std::get<7>(gradientC()), std::get<8>(
					vanLeer2()) * std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::vanLeer2Limiter::calculateNoRightLimit(
		const tensor3 & r, const tensor3 & gradientC) const noexcept
{
	tensor3 vanLeer2;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		vanLeer2.r()[j] = vanLeer2LimiterCalculation(r()[j]);

	return tensor3 { std::get<0>(vanLeer2()) * std::get<0>(gradientC()),
			std::get<1>(vanLeer2()) * std::get<1>(gradientC()), std::get<2>(
					vanLeer2()) * std::get<2>(gradientC()), std::get<3>(
					vanLeer2()) * std::get<3>(gradientC()), std::get<4>(
					vanLeer2()) * std::get<4>(gradientC()), std::get<5>(
					vanLeer2()) * std::get<5>(gradientC()), std::get<6>(
					vanLeer2()) * std::get<6>(gradientC()), std::get<7>(
					vanLeer2()) * std::get<7>(gradientC()), std::get<8>(
					vanLeer2()) * std::get<8>(gradientC()), std::get<9>(
					vanLeer2()) * std::get<9>(gradientC()), std::get<10>(
					vanLeer2()) * std::get<10>(gradientC()), std::get<11>(
					vanLeer2()) * std::get<11>(gradientC()), std::get<12>(
					vanLeer2()) * std::get<12>(gradientC()), std::get<13>(
					vanLeer2()) * std::get<13>(gradientC()), std::get<14>(
					vanLeer2()) * std::get<14>(gradientC()), std::get<15>(
					vanLeer2()) * std::get<15>(gradientC()), std::get<16>(
					vanLeer2()) * std::get<16>(gradientC()), std::get<17>(
					vanLeer2()) * std::get<17>(gradientC()), std::get<18>(
					vanLeer2()) * std::get<18>(gradientC()), std::get<19>(
					vanLeer2()) * std::get<19>(gradientC()), std::get<20>(
					vanLeer2()) * std::get<20>(gradientC()), std::get<21>(
					vanLeer2()) * std::get<21>(gradientC()), std::get<22>(
					vanLeer2()) * std::get<22>(gradientC()), std::get<23>(
					vanLeer2()) * std::get<23>(gradientC()), std::get<24>(
					vanLeer2()) * std::get<24>(gradientC()), std::get<25>(
					vanLeer2()) * std::get<25>(gradientC()), std::get<26>(
					vanLeer2()) * std::get<26>(gradientC()) };
}
