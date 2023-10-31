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
	vector HQUICK, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())) };

	for (std::size_t j = 0; j < vector::vsize; ++j)
		HQUICK.r()[j] = HQUICKLimiterCalculation(r()[j], xiR()[j]);

	return vector { std::get<0>(HQUICK()) * std::get<0>(gradientC()),
			std::get<1>(HQUICK()) * std::get<1>(gradientC()), std::get<2>(
					HQUICK()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::HQUICKLimiter::calculate(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor HQUICK, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())), 2 / (1 + std::get<3>(r())), 2
					/ (1 + std::get<4>(r())), 2 / (1 + std::get<5>(r())), 2
					/ (1 + std::get<6>(r())), 2 / (1 + std::get<7>(r())), 2
					/ (1 + std::get<8>(r())) };

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		HQUICK.r()[j] = HQUICKLimiterCalculation(r()[j], xiR()[j]);

	return tensor { std::get<0>(HQUICK()) * std::get<0>(gradientC()),
			std::get<1>(HQUICK()) * std::get<1>(gradientC()), std::get<2>(
					HQUICK()) * std::get<2>(gradientC()), std::get<3>(HQUICK())
					* std::get<3>(gradientC()), std::get<4>(HQUICK())
					* std::get<4>(gradientC()), std::get<5>(HQUICK())
					* std::get<5>(gradientC()), std::get<6>(HQUICK())
					* std::get<6>(gradientC()), std::get<7>(HQUICK())
					* std::get<7>(gradientC()), std::get<8>(HQUICK())
					* std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::HQUICKLimiter::calculate(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 HQUICK, xiR { 2 / (1 + std::get<0>(r())), 2
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
		HQUICK.r()[j] = HQUICKLimiterCalculation(r()[j], xiR()[j]);

	return tensor3 { std::get<0>(HQUICK()) * std::get<0>(gradientC()), std::get<
			1>(HQUICK()) * std::get<1>(gradientC()), std::get<2>(HQUICK())
			* std::get<2>(gradientC()), std::get<3>(HQUICK())
			* std::get<3>(gradientC()), std::get<4>(HQUICK())
			* std::get<4>(gradientC()), std::get<5>(HQUICK())
			* std::get<5>(gradientC()), std::get<6>(HQUICK())
			* std::get<6>(gradientC()), std::get<7>(HQUICK())
			* std::get<7>(gradientC()), std::get<8>(HQUICK())
			* std::get<8>(gradientC()), std::get<9>(HQUICK())
			* std::get<9>(gradientC()), std::get<10>(HQUICK())
			* std::get<10>(gradientC()), std::get<11>(HQUICK())
			* std::get<11>(gradientC()), std::get<12>(HQUICK())
			* std::get<12>(gradientC()), std::get<13>(HQUICK())
			* std::get<13>(gradientC()), std::get<14>(HQUICK())
			* std::get<14>(gradientC()), std::get<15>(HQUICK())
			* std::get<15>(gradientC()), std::get<16>(HQUICK())
			* std::get<16>(gradientC()), std::get<17>(HQUICK())
			* std::get<17>(gradientC()), std::get<18>(HQUICK())
			* std::get<18>(gradientC()), std::get<19>(HQUICK())
			* std::get<19>(gradientC()), std::get<20>(HQUICK())
			* std::get<20>(gradientC()), std::get<21>(HQUICK())
			* std::get<21>(gradientC()), std::get<22>(HQUICK())
			* std::get<22>(gradientC()), std::get<23>(HQUICK())
			* std::get<23>(gradientC()), std::get<24>(HQUICK())
			* std::get<24>(gradientC()), std::get<25>(HQUICK())
			* std::get<25>(gradientC()), std::get<26>(HQUICK())
			* std::get<26>(gradientC()) };
}

schemi::vector schemi::HQUICKLimiter::calculateNoRightLimit(const vector & r,
		const vector & gradientC) const noexcept
{
	vector HQUICK;

	for (std::size_t j = 0; j < vector::vsize; ++j)
		HQUICK.r()[j] = HQUICKLimiterCalculation(r()[j]);

	return vector { std::get<0>(HQUICK()) * std::get<0>(gradientC()),
			std::get<1>(HQUICK()) * std::get<1>(gradientC()), std::get<2>(
					HQUICK()) * std::get<2>(gradientC()) };
}

schemi::tensor schemi::HQUICKLimiter::calculateNoRightLimit(const tensor & r,
		const tensor & gradientC) const noexcept
{
	tensor HQUICK;

	for (std::size_t j = 0; j < tensor::vsize; ++j)
		HQUICK.r()[j] = HQUICKLimiterCalculation(r()[j]);

	return tensor { std::get<0>(HQUICK()) * std::get<0>(gradientC()),
			std::get<1>(HQUICK()) * std::get<1>(gradientC()), std::get<2>(
					HQUICK()) * std::get<2>(gradientC()), std::get<3>(HQUICK())
					* std::get<3>(gradientC()), std::get<4>(HQUICK())
					* std::get<4>(gradientC()), std::get<5>(HQUICK())
					* std::get<5>(gradientC()), std::get<6>(HQUICK())
					* std::get<6>(gradientC()), std::get<7>(HQUICK())
					* std::get<7>(gradientC()), std::get<8>(HQUICK())
					* std::get<8>(gradientC()) };
}

schemi::tensor3 schemi::HQUICKLimiter::calculateNoRightLimit(const tensor3 & r,
		const tensor3 & gradientC) const noexcept
{
	tensor3 HQUICK;

	for (std::size_t j = 0; j < tensor3::vsize; ++j)
		HQUICK.r()[j] = HQUICKLimiterCalculation(r()[j]);

	return tensor3 { std::get<0>(HQUICK()) * std::get<0>(gradientC()), std::get<
			1>(HQUICK()) * std::get<1>(gradientC()), std::get<2>(HQUICK())
			* std::get<2>(gradientC()), std::get<3>(HQUICK())
			* std::get<3>(gradientC()), std::get<4>(HQUICK())
			* std::get<4>(gradientC()), std::get<5>(HQUICK())
			* std::get<5>(gradientC()), std::get<6>(HQUICK())
			* std::get<6>(gradientC()), std::get<7>(HQUICK())
			* std::get<7>(gradientC()), std::get<8>(HQUICK())
			* std::get<8>(gradientC()), std::get<9>(HQUICK())
			* std::get<9>(gradientC()), std::get<10>(HQUICK())
			* std::get<10>(gradientC()), std::get<11>(HQUICK())
			* std::get<11>(gradientC()), std::get<12>(HQUICK())
			* std::get<12>(gradientC()), std::get<13>(HQUICK())
			* std::get<13>(gradientC()), std::get<14>(HQUICK())
			* std::get<14>(gradientC()), std::get<15>(HQUICK())
			* std::get<15>(gradientC()), std::get<16>(HQUICK())
			* std::get<16>(gradientC()), std::get<17>(HQUICK())
			* std::get<17>(gradientC()), std::get<18>(HQUICK())
			* std::get<18>(gradientC()), std::get<19>(HQUICK())
			* std::get<19>(gradientC()), std::get<20>(HQUICK())
			* std::get<20>(gradientC()), std::get<21>(HQUICK())
			* std::get<21>(gradientC()), std::get<22>(HQUICK())
			* std::get<22>(gradientC()), std::get<23>(HQUICK())
			* std::get<23>(gradientC()), std::get<24>(HQUICK())
			* std::get<24>(gradientC()), std::get<25>(HQUICK())
			* std::get<25>(gradientC()), std::get<26>(HQUICK())
			* std::get<26>(gradientC()) };
}
