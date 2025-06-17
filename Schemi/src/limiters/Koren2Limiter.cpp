/*
 * Koren2Limiter.cpp
 *
 *  Created on: 2025/05/23
 *      Author: Maxim Boldyrev
 */

#include "Koren2Limiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"

schemi::scalar schemi::Koren2Limiter::Koren2LimiterCalculation(const scalar r,
		const scalar xiR) const noexcept
{
	scalar Koren2;

	if (r <= 0)
		Koren2 = std::min(static_cast<scalar>(0),
				std::max(-2.0 / 3.0, (1 + 2 * r) / 3));
	else
		Koren2 = std::min(
				{ 4.0 * r / 3.0, std::min((1 + 2 * r) / 3,
						static_cast<scalar>(2)), 2 * xiR });

	return Koren2;
}

schemi::scalar schemi::Koren2Limiter::Koren2LimiterCalculation(
		const scalar r) const noexcept
{
	scalar Koren2;

	if (r <= 0)
		Koren2 = std::min(static_cast<scalar>(0),
				std::max(-2.0 / 3.0, (1 + 2 * r) / 3));
	else
		Koren2 = std::min(4.0 * r / 3.0,
				std::min((1 + 2 * r) / 3, static_cast<scalar>(2)));

	return Koren2;
}

schemi::vector schemi::Koren2Limiter::calculate(const vector & r,
		const vector & gradient) const noexcept
{
	vector Koren2, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), Koren2.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->Koren2LimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(Koren2()) * std::get<0>(gradient()),
			std::get<1>(Koren2()) * std::get<1>(gradient()), std::get<2>(
					Koren2()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::Koren2Limiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor Koren2, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())), 2 / (1 + std::get<3>(r())), 2
					/ (1 + std::get<4>(r())), 2 / (1 + std::get<5>(r())), 2
					/ (1 + std::get<6>(r())), 2 / (1 + std::get<7>(r())), 2
					/ (1 + std::get<8>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), Koren2.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->Koren2LimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(Koren2()) * std::get<0>(gradient()),
			std::get<1>(Koren2()) * std::get<1>(gradient()), std::get<2>(
					Koren2()) * std::get<2>(gradient()), std::get<3>(Koren2())
					* std::get<3>(gradient()), std::get<4>(Koren2())
					* std::get<4>(gradient()), std::get<5>(Koren2())
					* std::get<5>(gradient()), std::get<6>(Koren2())
					* std::get<6>(gradient()), std::get<7>(Koren2())
					* std::get<7>(gradient()), std::get<8>(Koren2())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::Koren2Limiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 Koren2, xiR { 2 / (1 + std::get<0>(r())), 2
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

	std::transform(r().begin(), r().end(), xiR().begin(), Koren2.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->Koren2LimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(Koren2()) * std::get<0>(gradient()),
			std::get<1>(Koren2()) * std::get<1>(gradient()), std::get<2>(
					Koren2()) * std::get<2>(gradient()), std::get<3>(Koren2())
					* std::get<3>(gradient()), std::get<4>(Koren2())
					* std::get<4>(gradient()), std::get<5>(Koren2())
					* std::get<5>(gradient()), std::get<6>(Koren2())
					* std::get<6>(gradient()), std::get<7>(Koren2())
					* std::get<7>(gradient()), std::get<8>(Koren2())
					* std::get<8>(gradient()), std::get<9>(Koren2())
					* std::get<9>(gradient()), std::get<10>(Koren2())
					* std::get<10>(gradient()), std::get<11>(Koren2())
					* std::get<11>(gradient()), std::get<12>(Koren2())
					* std::get<12>(gradient()), std::get<13>(Koren2())
					* std::get<13>(gradient()), std::get<14>(Koren2())
					* std::get<14>(gradient()), std::get<15>(Koren2())
					* std::get<15>(gradient()), std::get<16>(Koren2())
					* std::get<16>(gradient()), std::get<17>(Koren2())
					* std::get<17>(gradient()), std::get<18>(Koren2())
					* std::get<18>(gradient()), std::get<19>(Koren2())
					* std::get<19>(gradient()), std::get<20>(Koren2())
					* std::get<20>(gradient()), std::get<21>(Koren2())
					* std::get<21>(gradient()), std::get<22>(Koren2())
					* std::get<22>(gradient()), std::get<23>(Koren2())
					* std::get<23>(gradient()), std::get<24>(Koren2())
					* std::get<24>(gradient()), std::get<25>(Koren2())
					* std::get<25>(gradient()), std::get<26>(Koren2())
					* std::get<26>(gradient()) };
}

schemi::vector schemi::Koren2Limiter::calculateNoRSLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector Koren2;

	std::transform(r().begin(), r().end(), Koren2.r().begin(),
			[this](const auto r_j) 
			{	return this->Koren2LimiterCalculation(r_j);});

	return vector { std::get<0>(Koren2()) * std::get<0>(gradient()),
			std::get<1>(Koren2()) * std::get<1>(gradient()), std::get<2>(
					Koren2()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::Koren2Limiter::calculateNoRSLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor Koren2;

	std::transform(r().begin(), r().end(), Koren2.r().begin(),
			[this](const auto r_j) 
			{	return this->Koren2LimiterCalculation(r_j);});

	return tensor { std::get<0>(Koren2()) * std::get<0>(gradient()),
			std::get<1>(Koren2()) * std::get<1>(gradient()), std::get<2>(
					Koren2()) * std::get<2>(gradient()), std::get<3>(Koren2())
					* std::get<3>(gradient()), std::get<4>(Koren2())
					* std::get<4>(gradient()), std::get<5>(Koren2())
					* std::get<5>(gradient()), std::get<6>(Koren2())
					* std::get<6>(gradient()), std::get<7>(Koren2())
					* std::get<7>(gradient()), std::get<8>(Koren2())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::Koren2Limiter::calculateNoRSLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 Koren2;

	std::transform(r().begin(), r().end(), Koren2.r().begin(),
			[this](const auto r_j) 
			{	return this->Koren2LimiterCalculation(r_j);});

	return tensor3 { std::get<0>(Koren2()) * std::get<0>(gradient()),
			std::get<1>(Koren2()) * std::get<1>(gradient()), std::get<2>(
					Koren2()) * std::get<2>(gradient()), std::get<3>(Koren2())
					* std::get<3>(gradient()), std::get<4>(Koren2())
					* std::get<4>(gradient()), std::get<5>(Koren2())
					* std::get<5>(gradient()), std::get<6>(Koren2())
					* std::get<6>(gradient()), std::get<7>(Koren2())
					* std::get<7>(gradient()), std::get<8>(Koren2())
					* std::get<8>(gradient()), std::get<9>(Koren2())
					* std::get<9>(gradient()), std::get<10>(Koren2())
					* std::get<10>(gradient()), std::get<11>(Koren2())
					* std::get<11>(gradient()), std::get<12>(Koren2())
					* std::get<12>(gradient()), std::get<13>(Koren2())
					* std::get<13>(gradient()), std::get<14>(Koren2())
					* std::get<14>(gradient()), std::get<15>(Koren2())
					* std::get<15>(gradient()), std::get<16>(Koren2())
					* std::get<16>(gradient()), std::get<17>(Koren2())
					* std::get<17>(gradient()), std::get<18>(Koren2())
					* std::get<18>(gradient()), std::get<19>(Koren2())
					* std::get<19>(gradient()), std::get<20>(Koren2())
					* std::get<20>(gradient()), std::get<21>(Koren2())
					* std::get<21>(gradient()), std::get<22>(Koren2())
					* std::get<22>(gradient()), std::get<23>(Koren2())
					* std::get<23>(gradient()), std::get<24>(Koren2())
					* std::get<24>(gradient()), std::get<25>(Koren2())
					* std::get<25>(gradient()), std::get<26>(Koren2())
					* std::get<26>(gradient()) };
}
