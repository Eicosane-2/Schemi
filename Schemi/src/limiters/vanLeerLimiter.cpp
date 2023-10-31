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
	vector vanLeer, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanLeer.r()[j] = vanLeerLimiterCalculation(r()[j], xiR()[j]);

	return vector { std::get<0>(vanLeer()) * std::get<0>(gradientC()), std::get<
			1>(vanLeer()) * std::get<1>(gradientC()), std::get<2>(vanLeer())
			* std::get<2>(gradientC()) };
}

schemi::tensor schemi::vanLeerLimiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanLeer, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanLeer.r()[j] = vanLeerLimiterCalculation(r()[j], xiR()[j]);

	return tensor { std::get<0>(vanLeer()) * std::get<0>(gradientC()), std::get<
			1>(vanLeer()) * std::get<1>(gradientC()), std::get<2>(vanLeer())
			* std::get<2>(gradientC()), std::get<3>(vanLeer())
			* std::get<3>(gradientC()), std::get<4>(vanLeer())
			* std::get<4>(gradientC()), std::get<5>(vanLeer())
			* std::get<5>(gradientC()), std::get<6>(vanLeer())
			* std::get<6>(gradientC()), std::get<7>(vanLeer())
			* std::get<7>(gradientC()), std::get<8>(vanLeer())
			* std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::vanLeerLimiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 vanLeer, xiR { 2 / (1 + std::get<0>(r())), 2
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
		vanLeer.r()[j] = vanLeerLimiterCalculation(r()[j], xiR()[j]);

	return tensor3 { std::get<0>(vanLeer()) * std::get<0>(gradientC()),
			std::get<1>(vanLeer()) * std::get<1>(gradientC()), std::get<2>(
					vanLeer()) * std::get<2>(gradientC()), std::get<3>(
					vanLeer()) * std::get<3>(gradientC()), std::get<4>(
					vanLeer()) * std::get<4>(gradientC()), std::get<5>(
					vanLeer()) * std::get<5>(gradientC()), std::get<6>(
					vanLeer()) * std::get<6>(gradientC()), std::get<7>(
					vanLeer()) * std::get<7>(gradientC()), std::get<8>(
					vanLeer()) * std::get<8>(gradientC()), std::get<9>(
					vanLeer()) * std::get<9>(gradientC()), std::get<10>(
					vanLeer()) * std::get<10>(gradientC()), std::get<11>(
					vanLeer()) * std::get<11>(gradientC()), std::get<12>(
					vanLeer()) * std::get<12>(gradientC()), std::get<13>(
					vanLeer()) * std::get<13>(gradientC()), std::get<14>(
					vanLeer()) * std::get<14>(gradientC()), std::get<15>(
					vanLeer()) * std::get<15>(gradientC()), std::get<16>(
					vanLeer()) * std::get<16>(gradientC()), std::get<17>(
					vanLeer()) * std::get<17>(gradientC()), std::get<18>(
					vanLeer()) * std::get<18>(gradientC()), std::get<19>(
					vanLeer()) * std::get<19>(gradientC()), std::get<20>(
					vanLeer()) * std::get<20>(gradientC()), std::get<21>(
					vanLeer()) * std::get<21>(gradientC()), std::get<22>(
					vanLeer()) * std::get<22>(gradientC()), std::get<23>(
					vanLeer()) * std::get<23>(gradientC()), std::get<24>(
					vanLeer()) * std::get<24>(gradientC()), std::get<25>(
					vanLeer()) * std::get<25>(gradientC()), std::get<26>(
					vanLeer()) * std::get<26>(gradientC()) };
}

schemi::vector schemi::vanLeerLimiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector vanLeer;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		vanLeer.r()[j] = vanLeerLimiterCalculation(r()[j]);

	return vector { std::get<0>(vanLeer()) * std::get<0>(gradientC()), std::get<
			1>(vanLeer()) * std::get<1>(gradientC()), std::get<2>(vanLeer())
			* std::get<2>(gradientC()) };
}

schemi::tensor schemi::vanLeerLimiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor vanLeer;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		vanLeer.r()[j] = vanLeerLimiterCalculation(r()[j]);

	return tensor { std::get<0>(vanLeer()) * std::get<0>(gradientC()), std::get<
			1>(vanLeer()) * std::get<1>(gradientC()), std::get<2>(vanLeer())
			* std::get<2>(gradientC()), std::get<3>(vanLeer())
			* std::get<3>(gradientC()), std::get<4>(vanLeer())
			* std::get<4>(gradientC()), std::get<5>(vanLeer())
			* std::get<5>(gradientC()), std::get<6>(vanLeer())
			* std::get<6>(gradientC()), std::get<7>(vanLeer())
			* std::get<7>(gradientC()), std::get<8>(vanLeer())
			* std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::vanLeerLimiter::calculateNoRightLimit(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 vanLeer;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		vanLeer.r()[j] = vanLeerLimiterCalculation(r()[j]);

	return tensor3 { std::get<0>(vanLeer()) * std::get<0>(gradientC()),
			std::get<1>(vanLeer()) * std::get<1>(gradientC()), std::get<2>(
					vanLeer()) * std::get<2>(gradientC()), std::get<3>(
					vanLeer()) * std::get<3>(gradientC()), std::get<4>(
					vanLeer()) * std::get<4>(gradientC()), std::get<5>(
					vanLeer()) * std::get<5>(gradientC()), std::get<6>(
					vanLeer()) * std::get<6>(gradientC()), std::get<7>(
					vanLeer()) * std::get<7>(gradientC()), std::get<8>(
					vanLeer()) * std::get<8>(gradientC()), std::get<9>(
					vanLeer()) * std::get<9>(gradientC()), std::get<10>(
					vanLeer()) * std::get<10>(gradientC()), std::get<11>(
					vanLeer()) * std::get<11>(gradientC()), std::get<12>(
					vanLeer()) * std::get<12>(gradientC()), std::get<13>(
					vanLeer()) * std::get<13>(gradientC()), std::get<14>(
					vanLeer()) * std::get<14>(gradientC()), std::get<15>(
					vanLeer()) * std::get<15>(gradientC()), std::get<16>(
					vanLeer()) * std::get<16>(gradientC()), std::get<17>(
					vanLeer()) * std::get<17>(gradientC()), std::get<18>(
					vanLeer()) * std::get<18>(gradientC()), std::get<19>(
					vanLeer()) * std::get<19>(gradientC()), std::get<20>(
					vanLeer()) * std::get<20>(gradientC()), std::get<21>(
					vanLeer()) * std::get<21>(gradientC()), std::get<22>(
					vanLeer()) * std::get<22>(gradientC()), std::get<23>(
					vanLeer()) * std::get<23>(gradientC()), std::get<24>(
					vanLeer()) * std::get<24>(gradientC()), std::get<25>(
					vanLeer()) * std::get<25>(gradientC()), std::get<26>(
					vanLeer()) * std::get<26>(gradientC()) };
}
