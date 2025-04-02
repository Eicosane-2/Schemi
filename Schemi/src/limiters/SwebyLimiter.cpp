/*
 * SwebyLimiter.cpp
 *
 *  Created on: 2023/11/05
 *      Author: Maxim Boldyrev
 */

#include "SwebyLimiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"

schemi::scalar schemi::SwebyLimiter::SwebyLimiterCalculation(const scalar r,
		const scalar xiR) const noexcept
{
	scalar Sweby;

	if (r <= 0.)
		Sweby = 0;
	else
		Sweby = std::min(std::min(b, b * r), xiR);

	return Sweby;
}

schemi::scalar schemi::SwebyLimiter::SwebyLimiterCalculation(
		const scalar r) const noexcept
{
	scalar Sweby;

	if (r <= 0.)
		Sweby = 0;
	else
		Sweby = std::min(b, b * r);

	return Sweby;
}

schemi::vector schemi::SwebyLimiter::calculate(const vector & r,
		const vector & gradient) const noexcept
{
	vector Sweby, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), Sweby.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->SwebyLimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(Sweby()) * std::get<0>(gradient()), std::get<1>(
			Sweby()) * std::get<1>(gradient()), std::get<2>(Sweby())
			* std::get<2>(gradient()) };
}

schemi::tensor schemi::SwebyLimiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor Sweby, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())), 2 / (1 + std::get<3>(r())), 2
					/ (1 + std::get<4>(r())), 2 / (1 + std::get<5>(r())), 2
					/ (1 + std::get<6>(r())), 2 / (1 + std::get<7>(r())), 2
					/ (1 + std::get<8>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), Sweby.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->SwebyLimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(Sweby()) * std::get<0>(gradient()), std::get<1>(
			Sweby()) * std::get<1>(gradient()), std::get<2>(Sweby())
			* std::get<2>(gradient()), std::get<3>(Sweby())
			* std::get<3>(gradient()), std::get<4>(Sweby())
			* std::get<4>(gradient()), std::get<5>(Sweby())
			* std::get<5>(gradient()), std::get<6>(Sweby())
			* std::get<6>(gradient()), std::get<7>(Sweby())
			* std::get<7>(gradient()), std::get<8>(Sweby())
			* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::SwebyLimiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 Sweby, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())), 2 / (1 + std::get<3>(r())), 2
					/ (1 + std::get<4>(r())), 2 / (1 + std::get<5>(r())), 2
					/ (1 + std::get<6>(r())), 2 / (1 + std::get<7>(r())), 2
					/ (1 + std::get<8>(r())), 2 / (1 + std::get<9>(r())), 2
					/ (1 + std::get<10>(r())), 2 / (1 + std::get<11>(r())), 2
					/ (1 + std::get<12>(r())), 2 / (1 + std::get<13>(r())), 2
					/ (1 + std::get<14>(r())), 2 / (1 + std::get<15>(r())), 2
					/ (1 + std::get<16>(r())), 2 / (1 + std::get<17>(r())), 2
					/ (1 + std::get<18>(r())), 2 / (1 + std::get<19>(r())), 2
					/ (1 + std::get<20>(r())), 2 / (1 + std::get<21>(r())), 2
					/ (1 + std::get<22>(r())), 2 / (1 + std::get<23>(r())), 2
					/ (1 + std::get<24>(r())), 2 / (1 + std::get<25>(r())), 2
					/ (1 + std::get<26>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), Sweby.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->SwebyLimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(Sweby()) * std::get<0>(gradient()),
			std::get<1>(Sweby()) * std::get<1>(gradient()), std::get<2>(Sweby())
					* std::get<2>(gradient()), std::get<3>(Sweby())
					* std::get<3>(gradient()), std::get<4>(Sweby())
					* std::get<4>(gradient()), std::get<5>(Sweby())
					* std::get<5>(gradient()), std::get<6>(Sweby())
					* std::get<6>(gradient()), std::get<7>(Sweby())
					* std::get<7>(gradient()), std::get<8>(Sweby())
					* std::get<8>(gradient()), std::get<9>(Sweby())
					* std::get<9>(gradient()), std::get<10>(Sweby())
					* std::get<10>(gradient()), std::get<11>(Sweby())
					* std::get<11>(gradient()), std::get<12>(Sweby())
					* std::get<12>(gradient()), std::get<13>(Sweby())
					* std::get<13>(gradient()), std::get<14>(Sweby())
					* std::get<14>(gradient()), std::get<15>(Sweby())
					* std::get<15>(gradient()), std::get<16>(Sweby())
					* std::get<16>(gradient()), std::get<17>(Sweby())
					* std::get<17>(gradient()), std::get<18>(Sweby())
					* std::get<18>(gradient()), std::get<19>(Sweby())
					* std::get<19>(gradient()), std::get<20>(Sweby())
					* std::get<20>(gradient()), std::get<21>(Sweby())
					* std::get<21>(gradient()), std::get<22>(Sweby())
					* std::get<22>(gradient()), std::get<23>(Sweby())
					* std::get<23>(gradient()), std::get<24>(Sweby())
					* std::get<24>(gradient()), std::get<25>(Sweby())
					* std::get<25>(gradient()), std::get<26>(Sweby())
					* std::get<26>(gradient()) };
}

schemi::vector schemi::SwebyLimiter::calculateNoRSLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector Sweby;

	std::transform(r().begin(), r().end(), Sweby.r().begin(),
			[this](const auto r_j) 
			{	return this->SwebyLimiterCalculation(r_j);});

	const auto beta = elementsDivision((vector(1) + 2 * r), vector(3));

	std::transform(Sweby().begin(), Sweby().end(), beta().begin(),
			Sweby.r().begin(), [](const auto limiter_j, const auto beta_j) 
			{	return std::max(std::min(limiter_j, beta_j),0.0);});

	return vector { std::get<0>(Sweby()) * std::get<0>(gradient()), std::get<1>(
			Sweby()) * std::get<1>(gradient()), std::get<2>(Sweby())
			* std::get<2>(gradient()) };
}

schemi::tensor schemi::SwebyLimiter::calculateNoRSLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor Sweby;

	std::transform(r().begin(), r().end(), Sweby.r().begin(),
			[this](const auto r_j) 
			{	return this->SwebyLimiterCalculation(r_j);});

	const auto beta = elementsDivision((tensor(1) + 2 * r), tensor(3));

	std::transform(Sweby().begin(), Sweby().end(), beta().begin(),
			Sweby.r().begin(), [](const auto limiter_j, const auto beta_j) 
			{	return std::max(std::min(limiter_j, beta_j),0.0);});

	return tensor { std::get<0>(Sweby()) * std::get<0>(gradient()), std::get<1>(
			Sweby()) * std::get<1>(gradient()), std::get<2>(Sweby())
			* std::get<2>(gradient()), std::get<3>(Sweby())
			* std::get<3>(gradient()), std::get<4>(Sweby())
			* std::get<4>(gradient()), std::get<5>(Sweby())
			* std::get<5>(gradient()), std::get<6>(Sweby())
			* std::get<6>(gradient()), std::get<7>(Sweby())
			* std::get<7>(gradient()), std::get<8>(Sweby())
			* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::SwebyLimiter::calculateNoRSLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 Sweby;

	std::transform(r().begin(), r().end(), Sweby.r().begin(),
			[this](const auto r_j) 
			{	return this->SwebyLimiterCalculation(r_j);});

	const auto beta = elementsDivision((tensor3(1) + 2 * r), tensor3(3));

	std::transform(Sweby().begin(), Sweby().end(), beta().begin(),
			Sweby.r().begin(), [](const auto limiter_j, const auto beta_j) 
			{	return std::max(std::min(limiter_j, beta_j),0.0);});

	return tensor3 { std::get<0>(Sweby()) * std::get<0>(gradient()),
			std::get<1>(Sweby()) * std::get<1>(gradient()), std::get<2>(Sweby())
					* std::get<2>(gradient()), std::get<3>(Sweby())
					* std::get<3>(gradient()), std::get<4>(Sweby())
					* std::get<4>(gradient()), std::get<5>(Sweby())
					* std::get<5>(gradient()), std::get<6>(Sweby())
					* std::get<6>(gradient()), std::get<7>(Sweby())
					* std::get<7>(gradient()), std::get<8>(Sweby())
					* std::get<8>(gradient()), std::get<9>(Sweby())
					* std::get<9>(gradient()), std::get<10>(Sweby())
					* std::get<10>(gradient()), std::get<11>(Sweby())
					* std::get<11>(gradient()), std::get<12>(Sweby())
					* std::get<12>(gradient()), std::get<13>(Sweby())
					* std::get<13>(gradient()), std::get<14>(Sweby())
					* std::get<14>(gradient()), std::get<15>(Sweby())
					* std::get<15>(gradient()), std::get<16>(Sweby())
					* std::get<16>(gradient()), std::get<17>(Sweby())
					* std::get<17>(gradient()), std::get<18>(Sweby())
					* std::get<18>(gradient()), std::get<19>(Sweby())
					* std::get<19>(gradient()), std::get<20>(Sweby())
					* std::get<20>(gradient()), std::get<21>(Sweby())
					* std::get<21>(gradient()), std::get<22>(Sweby())
					* std::get<22>(gradient()), std::get<23>(Sweby())
					* std::get<23>(gradient()), std::get<24>(Sweby())
					* std::get<24>(gradient()), std::get<25>(Sweby())
					* std::get<25>(gradient()), std::get<26>(Sweby())
					* std::get<26>(gradient()) };
}
