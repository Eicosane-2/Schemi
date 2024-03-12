/*
 * minmodLimiter.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "minmodLimiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"

schemi::scalar schemi::minmodLimiter::minmodLimiterCalculation(const scalar r,
		const scalar xiR) const noexcept
{
	scalar minmod;

	if (r <= 0.)
		minmod = 0;
	else if (r <= 1.)
		minmod = r;
	else
		minmod = std::min(static_cast<scalar>(1.), xiR);

	return minmod;
}

schemi::scalar schemi::minmodLimiter::minmodLimiterCalculation(
		const scalar r) const noexcept
{
	scalar minmod;

	if (r <= 0.)
		minmod = 0;
	else if (r <= 1.)
		minmod = r;
	else
		minmod = 1.;

	return minmod;
}

schemi::vector schemi::minmodLimiter::calculate(const vector & r,
		const vector & gradient) const noexcept
{
	vector minmod, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), minmod.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->minmodLimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(minmod()) * std::get<0>(gradient()),
			std::get<1>(minmod()) * std::get<1>(gradient()), std::get<2>(
					minmod()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::minmodLimiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor minmod, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())), 2 / (1 + std::get<3>(r())), 2
					/ (1 + std::get<4>(r())), 2 / (1 + std::get<5>(r())), 2
					/ (1 + std::get<6>(r())), 2 / (1 + std::get<7>(r())), 2
					/ (1 + std::get<8>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), minmod.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->minmodLimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(minmod()) * std::get<0>(gradient()),
			std::get<1>(minmod()) * std::get<1>(gradient()), std::get<2>(
					minmod()) * std::get<2>(gradient()), std::get<3>(minmod())
					* std::get<3>(gradient()), std::get<4>(minmod())
					* std::get<4>(gradient()), std::get<5>(minmod())
					* std::get<5>(gradient()), std::get<6>(minmod())
					* std::get<6>(gradient()), std::get<7>(minmod())
					* std::get<7>(gradient()), std::get<8>(minmod())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::minmodLimiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 minmod, xiR { 2 / (1 + std::get<0>(r())), 2
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

	std::transform(r().begin(), r().end(), xiR().begin(), minmod.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->minmodLimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(minmod()) * std::get<0>(gradient()),
			std::get<1>(minmod()) * std::get<1>(gradient()), std::get<2>(
					minmod()) * std::get<2>(gradient()), std::get<3>(minmod())
					* std::get<3>(gradient()), std::get<4>(minmod())
					* std::get<4>(gradient()), std::get<5>(minmod())
					* std::get<5>(gradient()), std::get<6>(minmod())
					* std::get<6>(gradient()), std::get<7>(minmod())
					* std::get<7>(gradient()), std::get<8>(minmod())
					* std::get<8>(gradient()), std::get<9>(minmod())
					* std::get<9>(gradient()), std::get<10>(minmod())
					* std::get<10>(gradient()), std::get<11>(minmod())
					* std::get<11>(gradient()), std::get<12>(minmod())
					* std::get<12>(gradient()), std::get<13>(minmod())
					* std::get<13>(gradient()), std::get<14>(minmod())
					* std::get<14>(gradient()), std::get<15>(minmod())
					* std::get<15>(gradient()), std::get<16>(minmod())
					* std::get<16>(gradient()), std::get<17>(minmod())
					* std::get<17>(gradient()), std::get<18>(minmod())
					* std::get<18>(gradient()), std::get<19>(minmod())
					* std::get<19>(gradient()), std::get<20>(minmod())
					* std::get<20>(gradient()), std::get<21>(minmod())
					* std::get<21>(gradient()), std::get<22>(minmod())
					* std::get<22>(gradient()), std::get<23>(minmod())
					* std::get<23>(gradient()), std::get<24>(minmod())
					* std::get<24>(gradient()), std::get<25>(minmod())
					* std::get<25>(gradient()), std::get<26>(minmod())
					* std::get<26>(gradient()) };
}

schemi::vector schemi::minmodLimiter::calculate3OLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector minmod;

	std::transform(r().begin(), r().end(), minmod.r().begin(),
			[this](const auto r_j) 
			{	return this->minmodLimiterCalculation(r_j);});

	const auto beta = elementsDivision((vector(1) + 2 * r), vector(3));

	std::transform(minmod().begin(), minmod().end(), beta().begin(),
			minmod.r().begin(), [](const auto limiter_j, const auto beta_j) 
			{	return std::max(std::min(limiter_j, beta_j),0.0);});

	return vector { std::get<0>(minmod()) * std::get<0>(gradient()),
			std::get<1>(minmod()) * std::get<1>(gradient()), std::get<2>(
					minmod()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::minmodLimiter::calculate3OLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor minmod;

	std::transform(r().begin(), r().end(), minmod.r().begin(),
			[this](const auto r_j) 
			{	return this->minmodLimiterCalculation(r_j);});

	const auto beta = elementsDivision((tensor(1) + 2 * r), tensor(3));

	std::transform(minmod().begin(), minmod().end(), beta().begin(),
			minmod.r().begin(), [](const auto limiter_j, const auto beta_j) 
			{	return std::max(std::min(limiter_j, beta_j),0.0);});

	return tensor { std::get<0>(minmod()) * std::get<0>(gradient()),
			std::get<1>(minmod()) * std::get<1>(gradient()), std::get<2>(
					minmod()) * std::get<2>(gradient()), std::get<3>(minmod())
					* std::get<3>(gradient()), std::get<4>(minmod())
					* std::get<4>(gradient()), std::get<5>(minmod())
					* std::get<5>(gradient()), std::get<6>(minmod())
					* std::get<6>(gradient()), std::get<7>(minmod())
					* std::get<7>(gradient()), std::get<8>(minmod())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::minmodLimiter::calculate3OLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 minmod;

	std::transform(r().begin(), r().end(), minmod.r().begin(),
			[this](const auto r_j) 
			{	return this->minmodLimiterCalculation(r_j);});

	const auto beta = elementsDivision((tensor3(1) + 2 * r), tensor3(3));

	std::transform(minmod().begin(), minmod().end(), beta().begin(),
			minmod.r().begin(), [](const auto limiter_j, const auto beta_j) 
			{	return std::max(std::min(limiter_j, beta_j),0.0);});

	return tensor3 { std::get<0>(minmod()) * std::get<0>(gradient()),
			std::get<1>(minmod()) * std::get<1>(gradient()), std::get<2>(
					minmod()) * std::get<2>(gradient()), std::get<3>(minmod())
					* std::get<3>(gradient()), std::get<4>(minmod())
					* std::get<4>(gradient()), std::get<5>(minmod())
					* std::get<5>(gradient()), std::get<6>(minmod())
					* std::get<6>(gradient()), std::get<7>(minmod())
					* std::get<7>(gradient()), std::get<8>(minmod())
					* std::get<8>(gradient()), std::get<9>(minmod())
					* std::get<9>(gradient()), std::get<10>(minmod())
					* std::get<10>(gradient()), std::get<11>(minmod())
					* std::get<11>(gradient()), std::get<12>(minmod())
					* std::get<12>(gradient()), std::get<13>(minmod())
					* std::get<13>(gradient()), std::get<14>(minmod())
					* std::get<14>(gradient()), std::get<15>(minmod())
					* std::get<15>(gradient()), std::get<16>(minmod())
					* std::get<16>(gradient()), std::get<17>(minmod())
					* std::get<17>(gradient()), std::get<18>(minmod())
					* std::get<18>(gradient()), std::get<19>(minmod())
					* std::get<19>(gradient()), std::get<20>(minmod())
					* std::get<20>(gradient()), std::get<21>(minmod())
					* std::get<21>(gradient()), std::get<22>(minmod())
					* std::get<22>(gradient()), std::get<23>(minmod())
					* std::get<23>(gradient()), std::get<24>(minmod())
					* std::get<24>(gradient()), std::get<25>(minmod())
					* std::get<25>(gradient()), std::get<26>(minmod())
					* std::get<26>(gradient()) };
}
