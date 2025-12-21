/*
 * HQUICKLimiter.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HQUICKLimiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"

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
		const vector & gradient) const noexcept
{
	vector HQUICK, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())) };

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			HQUICK.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->HQUICKLimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(HQUICK()) * std::get<0>(gradient()),
			std::get<1>(HQUICK()) * std::get<1>(gradient()), std::get<2>(
					HQUICK()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::HQUICKLimiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor HQUICK, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())), 2 / (1 + std::get<3>(r())), 2
					/ (1 + std::get<4>(r())), 2 / (1 + std::get<5>(r())), 2
					/ (1 + std::get<6>(r())), 2 / (1 + std::get<7>(r())), 2
					/ (1 + std::get<8>(r())) };

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			HQUICK.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->HQUICKLimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(HQUICK()) * std::get<0>(gradient()),
			std::get<1>(HQUICK()) * std::get<1>(gradient()), std::get<2>(
					HQUICK()) * std::get<2>(gradient()), std::get<3>(HQUICK())
					* std::get<3>(gradient()), std::get<4>(HQUICK())
					* std::get<4>(gradient()), std::get<5>(HQUICK())
					* std::get<5>(gradient()), std::get<6>(HQUICK())
					* std::get<6>(gradient()), std::get<7>(HQUICK())
					* std::get<7>(gradient()), std::get<8>(HQUICK())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::HQUICKLimiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
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

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			HQUICK.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->HQUICKLimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(HQUICK()) * std::get<0>(gradient()),
			std::get<1>(HQUICK()) * std::get<1>(gradient()), std::get<2>(
					HQUICK()) * std::get<2>(gradient()), std::get<3>(HQUICK())
					* std::get<3>(gradient()), std::get<4>(HQUICK())
					* std::get<4>(gradient()), std::get<5>(HQUICK())
					* std::get<5>(gradient()), std::get<6>(HQUICK())
					* std::get<6>(gradient()), std::get<7>(HQUICK())
					* std::get<7>(gradient()), std::get<8>(HQUICK())
					* std::get<8>(gradient()), std::get<9>(HQUICK())
					* std::get<9>(gradient()), std::get<10>(HQUICK())
					* std::get<10>(gradient()), std::get<11>(HQUICK())
					* std::get<11>(gradient()), std::get<12>(HQUICK())
					* std::get<12>(gradient()), std::get<13>(HQUICK())
					* std::get<13>(gradient()), std::get<14>(HQUICK())
					* std::get<14>(gradient()), std::get<15>(HQUICK())
					* std::get<15>(gradient()), std::get<16>(HQUICK())
					* std::get<16>(gradient()), std::get<17>(HQUICK())
					* std::get<17>(gradient()), std::get<18>(HQUICK())
					* std::get<18>(gradient()), std::get<19>(HQUICK())
					* std::get<19>(gradient()), std::get<20>(HQUICK())
					* std::get<20>(gradient()), std::get<21>(HQUICK())
					* std::get<21>(gradient()), std::get<22>(HQUICK())
					* std::get<22>(gradient()), std::get<23>(HQUICK())
					* std::get<23>(gradient()), std::get<24>(HQUICK())
					* std::get<24>(gradient()), std::get<25>(HQUICK())
					* std::get<25>(gradient()), std::get<26>(HQUICK())
					* std::get<26>(gradient()) };
}

schemi::vector schemi::HQUICKLimiter::calculateNoRSLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector HQUICK;

	std::transform(r().cbegin(), r().cend(), HQUICK.wr().begin(),
			[this](const auto r_j) 
			{	return this->HQUICKLimiterCalculation(r_j);});

	return vector { std::get<0>(HQUICK()) * std::get<0>(gradient()),
			std::get<1>(HQUICK()) * std::get<1>(gradient()), std::get<2>(
					HQUICK()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::HQUICKLimiter::calculateNoRSLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor HQUICK;

	std::transform(r().cbegin(), r().cend(), HQUICK.wr().begin(),
			[this](const auto r_j) 
			{	return this->HQUICKLimiterCalculation(r_j);});

	return tensor { std::get<0>(HQUICK()) * std::get<0>(gradient()),
			std::get<1>(HQUICK()) * std::get<1>(gradient()), std::get<2>(
					HQUICK()) * std::get<2>(gradient()), std::get<3>(HQUICK())
					* std::get<3>(gradient()), std::get<4>(HQUICK())
					* std::get<4>(gradient()), std::get<5>(HQUICK())
					* std::get<5>(gradient()), std::get<6>(HQUICK())
					* std::get<6>(gradient()), std::get<7>(HQUICK())
					* std::get<7>(gradient()), std::get<8>(HQUICK())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::HQUICKLimiter::calculateNoRSLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 HQUICK;

	std::transform(r().cbegin(), r().cend(), HQUICK.wr().begin(),
			[this](const auto r_j) 
			{	return this->HQUICKLimiterCalculation(r_j);});

	return tensor3 { std::get<0>(HQUICK()) * std::get<0>(gradient()),
			std::get<1>(HQUICK()) * std::get<1>(gradient()), std::get<2>(
					HQUICK()) * std::get<2>(gradient()), std::get<3>(HQUICK())
					* std::get<3>(gradient()), std::get<4>(HQUICK())
					* std::get<4>(gradient()), std::get<5>(HQUICK())
					* std::get<5>(gradient()), std::get<6>(HQUICK())
					* std::get<6>(gradient()), std::get<7>(HQUICK())
					* std::get<7>(gradient()), std::get<8>(HQUICK())
					* std::get<8>(gradient()), std::get<9>(HQUICK())
					* std::get<9>(gradient()), std::get<10>(HQUICK())
					* std::get<10>(gradient()), std::get<11>(HQUICK())
					* std::get<11>(gradient()), std::get<12>(HQUICK())
					* std::get<12>(gradient()), std::get<13>(HQUICK())
					* std::get<13>(gradient()), std::get<14>(HQUICK())
					* std::get<14>(gradient()), std::get<15>(HQUICK())
					* std::get<15>(gradient()), std::get<16>(HQUICK())
					* std::get<16>(gradient()), std::get<17>(HQUICK())
					* std::get<17>(gradient()), std::get<18>(HQUICK())
					* std::get<18>(gradient()), std::get<19>(HQUICK())
					* std::get<19>(gradient()), std::get<20>(HQUICK())
					* std::get<20>(gradient()), std::get<21>(HQUICK())
					* std::get<21>(gradient()), std::get<22>(HQUICK())
					* std::get<22>(gradient()), std::get<23>(HQUICK())
					* std::get<23>(gradient()), std::get<24>(HQUICK())
					* std::get<24>(gradient()), std::get<25>(HQUICK())
					* std::get<25>(gradient()), std::get<26>(HQUICK())
					* std::get<26>(gradient()) };
}
