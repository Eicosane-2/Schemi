/*
 * vanAlbada2Limiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vanAlbada2Limiter.hpp"

#include "intExpPow.hpp"

schemi::scalar schemi::vanAlbada2Limiter::vanAlbada2LimiterCalculation(
		const scalar r, const scalar xiR) const noexcept
{
	scalar vanAlbada2;

	if (r <= 0.)
		vanAlbada2 = 0;
	else
		vanAlbada2 = std::min(2 * r / (1 + pow<scalar, 2>(r)), 2 * xiR);

	return vanAlbada2;
}

schemi::scalar schemi::vanAlbada2Limiter::vanAlbada2LimiterCalculation(
		const scalar r) const noexcept
{
	scalar vanAlbada2;

	if (r <= 0.)
		vanAlbada2 = 0;
	else
		vanAlbada2 = 2 * r / (1 + pow<scalar, 2>(r));

	return vanAlbada2;
}

schemi::vector schemi::vanAlbada2Limiter::calculate(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanAlbada2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanAlbada2.r()[j] = vanAlbada2LimiterCalculation(r()[j], xiR()[j]);

	return vector { std::get<0>(vanAlbada2()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::vanAlbada2Limiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanAlbada2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanAlbada2.r()[j] = vanAlbada2LimiterCalculation(r()[j], xiR()[j]);

	return tensor { std::get<0>(vanAlbada2()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradientC()), std::get<3>(
					vanAlbada2()) * std::get<3>(gradientC()), std::get<4>(
					vanAlbada2()) * std::get<4>(gradientC()), std::get<5>(
					vanAlbada2()) * std::get<5>(gradientC()), std::get<6>(
					vanAlbada2()) * std::get<6>(gradientC()), std::get<7>(
					vanAlbada2()) * std::get<7>(gradientC()), std::get<8>(
					vanAlbada2()) * std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::vanAlbada2Limiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 vanAlbada2, xiR { 2 / (1 + std::get<0>(r())), 2
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
		vanAlbada2.r()[j] = vanAlbada2LimiterCalculation(r()[j], xiR()[j]);

	return tensor3 { std::get<0>(vanAlbada2()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradientC()), std::get<3>(
					vanAlbada2()) * std::get<3>(gradientC()), std::get<4>(
					vanAlbada2()) * std::get<4>(gradientC()), std::get<5>(
					vanAlbada2()) * std::get<5>(gradientC()), std::get<6>(
					vanAlbada2()) * std::get<6>(gradientC()), std::get<7>(
					vanAlbada2()) * std::get<7>(gradientC()), std::get<8>(
					vanAlbada2()) * std::get<8>(gradientC()), std::get<9>(
					vanAlbada2()) * std::get<9>(gradientC()), std::get<10>(
					vanAlbada2()) * std::get<10>(gradientC()), std::get<11>(
					vanAlbada2()) * std::get<11>(gradientC()), std::get<12>(
					vanAlbada2()) * std::get<12>(gradientC()), std::get<13>(
					vanAlbada2()) * std::get<13>(gradientC()), std::get<14>(
					vanAlbada2()) * std::get<14>(gradientC()), std::get<15>(
					vanAlbada2()) * std::get<15>(gradientC()), std::get<16>(
					vanAlbada2()) * std::get<16>(gradientC()), std::get<17>(
					vanAlbada2()) * std::get<17>(gradientC()), std::get<18>(
					vanAlbada2()) * std::get<18>(gradientC()), std::get<19>(
					vanAlbada2()) * std::get<19>(gradientC()), std::get<20>(
					vanAlbada2()) * std::get<20>(gradientC()), std::get<21>(
					vanAlbada2()) * std::get<21>(gradientC()), std::get<22>(
					vanAlbada2()) * std::get<22>(gradientC()), std::get<23>(
					vanAlbada2()) * std::get<23>(gradientC()), std::get<24>(
					vanAlbada2()) * std::get<24>(gradientC()), std::get<25>(
					vanAlbada2()) * std::get<25>(gradientC()), std::get<26>(
					vanAlbada2()) * std::get<26>(gradientC()) };
}

schemi::vector schemi::vanAlbada2Limiter::calculateNoRightLimit(
		const vector & r, const vector & gradientC) const noexcept
{
	vector vanAlbada2;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanAlbada2.r()[j] = vanAlbada2LimiterCalculation(r()[j]);

	return vector { std::get<0>(vanAlbada2()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::vanAlbada2Limiter::calculateNoRightLimit(
		const tensor & r, const tensor & gradientC) const noexcept
{
	tensor vanAlbada2;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanAlbada2.r()[j] = vanAlbada2LimiterCalculation(r()[j]);

	return tensor { std::get<0>(vanAlbada2()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradientC()), std::get<3>(
					vanAlbada2()) * std::get<3>(gradientC()), std::get<4>(
					vanAlbada2()) * std::get<4>(gradientC()), std::get<5>(
					vanAlbada2()) * std::get<5>(gradientC()), std::get<6>(
					vanAlbada2()) * std::get<6>(gradientC()), std::get<7>(
					vanAlbada2()) * std::get<7>(gradientC()), std::get<8>(
					vanAlbada2()) * std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::vanAlbada2Limiter::calculateNoRightLimit(
		const tensor3 & r, const tensor3 & gradientC) const noexcept
{
	tensor3 vanAlbada2;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		vanAlbada2.r()[j] = vanAlbada2LimiterCalculation(r()[j]);

	return tensor3 { std::get<0>(vanAlbada2()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradientC()), std::get<3>(
					vanAlbada2()) * std::get<3>(gradientC()), std::get<4>(
					vanAlbada2()) * std::get<4>(gradientC()), std::get<5>(
					vanAlbada2()) * std::get<5>(gradientC()), std::get<6>(
					vanAlbada2()) * std::get<6>(gradientC()), std::get<7>(
					vanAlbada2()) * std::get<7>(gradientC()), std::get<8>(
					vanAlbada2()) * std::get<8>(gradientC()), std::get<9>(
					vanAlbada2()) * std::get<9>(gradientC()), std::get<10>(
					vanAlbada2()) * std::get<10>(gradientC()), std::get<11>(
					vanAlbada2()) * std::get<11>(gradientC()), std::get<12>(
					vanAlbada2()) * std::get<12>(gradientC()), std::get<13>(
					vanAlbada2()) * std::get<13>(gradientC()), std::get<14>(
					vanAlbada2()) * std::get<14>(gradientC()), std::get<15>(
					vanAlbada2()) * std::get<15>(gradientC()), std::get<16>(
					vanAlbada2()) * std::get<16>(gradientC()), std::get<17>(
					vanAlbada2()) * std::get<17>(gradientC()), std::get<18>(
					vanAlbada2()) * std::get<18>(gradientC()), std::get<19>(
					vanAlbada2()) * std::get<19>(gradientC()), std::get<20>(
					vanAlbada2()) * std::get<20>(gradientC()), std::get<21>(
					vanAlbada2()) * std::get<21>(gradientC()), std::get<22>(
					vanAlbada2()) * std::get<22>(gradientC()), std::get<23>(
					vanAlbada2()) * std::get<23>(gradientC()), std::get<24>(
					vanAlbada2()) * std::get<24>(gradientC()), std::get<25>(
					vanAlbada2()) * std::get<25>(gradientC()), std::get<26>(
					vanAlbada2()) * std::get<26>(gradientC()) };
}
