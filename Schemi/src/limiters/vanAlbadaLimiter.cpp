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
	vector vanAlbada, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanAlbada.r()[j] = vanAlbadaLimiterCalculation(r()[j], xiR()[j]);

	return vector { std::get<0>(vanAlbada()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::vanAlbadaLimiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanAlbada, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanAlbada.r()[j] = vanAlbadaLimiterCalculation(r()[j], xiR()[j]);

	return tensor { std::get<0>(vanAlbada()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada()) * std::get<2>(gradientC()), std::get<3>(
					vanAlbada()) * std::get<3>(gradientC()), std::get<4>(
					vanAlbada()) * std::get<4>(gradientC()), std::get<5>(
					vanAlbada()) * std::get<5>(gradientC()), std::get<6>(
					vanAlbada()) * std::get<6>(gradientC()), std::get<7>(
					vanAlbada()) * std::get<7>(gradientC()), std::get<8>(
					vanAlbada()) * std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::vanAlbadaLimiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 vanAlbada, xiR { 2 / (1 + std::get<0>(r())), 2
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
		vanAlbada.r()[j] = vanAlbadaLimiterCalculation(r()[j], xiR()[j]);

	return tensor3 { std::get<0>(vanAlbada()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada()) * std::get<2>(gradientC()), std::get<3>(
					vanAlbada()) * std::get<3>(gradientC()), std::get<4>(
					vanAlbada()) * std::get<4>(gradientC()), std::get<5>(
					vanAlbada()) * std::get<5>(gradientC()), std::get<6>(
					vanAlbada()) * std::get<6>(gradientC()), std::get<7>(
					vanAlbada()) * std::get<7>(gradientC()), std::get<8>(
					vanAlbada()) * std::get<8>(gradientC()), std::get<9>(
					vanAlbada()) * std::get<9>(gradientC()), std::get<10>(
					vanAlbada()) * std::get<10>(gradientC()), std::get<11>(
					vanAlbada()) * std::get<11>(gradientC()), std::get<12>(
					vanAlbada()) * std::get<12>(gradientC()), std::get<13>(
					vanAlbada()) * std::get<13>(gradientC()), std::get<14>(
					vanAlbada()) * std::get<14>(gradientC()), std::get<15>(
					vanAlbada()) * std::get<15>(gradientC()), std::get<16>(
					vanAlbada()) * std::get<16>(gradientC()), std::get<17>(
					vanAlbada()) * std::get<17>(gradientC()), std::get<18>(
					vanAlbada()) * std::get<18>(gradientC()), std::get<19>(
					vanAlbada()) * std::get<19>(gradientC()), std::get<20>(
					vanAlbada()) * std::get<20>(gradientC()), std::get<21>(
					vanAlbada()) * std::get<21>(gradientC()), std::get<22>(
					vanAlbada()) * std::get<22>(gradientC()), std::get<23>(
					vanAlbada()) * std::get<23>(gradientC()), std::get<24>(
					vanAlbada()) * std::get<24>(gradientC()), std::get<25>(
					vanAlbada()) * std::get<25>(gradientC()), std::get<26>(
					vanAlbada()) * std::get<26>(gradientC()) };
}

schemi::vector schemi::vanAlbadaLimiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanAlbada;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanAlbada.r()[j] = vanAlbadaLimiterCalculation(r()[j]);

	return vector { std::get<0>(vanAlbada()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::vanAlbadaLimiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanAlbada;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanAlbada.r()[j] = vanAlbadaLimiterCalculation(r()[j]);

	return tensor { std::get<0>(vanAlbada()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada()) * std::get<2>(gradientC()), std::get<3>(
					vanAlbada()) * std::get<3>(gradientC()), std::get<4>(
					vanAlbada()) * std::get<4>(gradientC()), std::get<5>(
					vanAlbada()) * std::get<5>(gradientC()), std::get<6>(
					vanAlbada()) * std::get<6>(gradientC()), std::get<7>(
					vanAlbada()) * std::get<7>(gradientC()), std::get<8>(
					vanAlbada()) * std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::vanAlbadaLimiter::calculateNoRightLimit(
		const tensor3 & r, const tensor3 & gradientC) const noexcept
{
	tensor3 vanAlbada;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		vanAlbada.r()[j] = vanAlbadaLimiterCalculation(r()[j]);

	return tensor3 { std::get<0>(vanAlbada()) * std::get<0>(gradientC()),
			std::get<1>(vanAlbada()) * std::get<1>(gradientC()), std::get<2>(
					vanAlbada()) * std::get<2>(gradientC()), std::get<3>(
					vanAlbada()) * std::get<3>(gradientC()), std::get<4>(
					vanAlbada()) * std::get<4>(gradientC()), std::get<5>(
					vanAlbada()) * std::get<5>(gradientC()), std::get<6>(
					vanAlbada()) * std::get<6>(gradientC()), std::get<7>(
					vanAlbada()) * std::get<7>(gradientC()), std::get<8>(
					vanAlbada()) * std::get<8>(gradientC()), std::get<9>(
					vanAlbada()) * std::get<9>(gradientC()), std::get<10>(
					vanAlbada()) * std::get<10>(gradientC()), std::get<11>(
					vanAlbada()) * std::get<11>(gradientC()), std::get<12>(
					vanAlbada()) * std::get<12>(gradientC()), std::get<13>(
					vanAlbada()) * std::get<13>(gradientC()), std::get<14>(
					vanAlbada()) * std::get<14>(gradientC()), std::get<15>(
					vanAlbada()) * std::get<15>(gradientC()), std::get<16>(
					vanAlbada()) * std::get<16>(gradientC()), std::get<17>(
					vanAlbada()) * std::get<17>(gradientC()), std::get<18>(
					vanAlbada()) * std::get<18>(gradientC()), std::get<19>(
					vanAlbada()) * std::get<19>(gradientC()), std::get<20>(
					vanAlbada()) * std::get<20>(gradientC()), std::get<21>(
					vanAlbada()) * std::get<21>(gradientC()), std::get<22>(
					vanAlbada()) * std::get<22>(gradientC()), std::get<23>(
					vanAlbada()) * std::get<23>(gradientC()), std::get<24>(
					vanAlbada()) * std::get<24>(gradientC()), std::get<25>(
					vanAlbada()) * std::get<25>(gradientC()), std::get<26>(
					vanAlbada()) * std::get<26>(gradientC()) };
}
