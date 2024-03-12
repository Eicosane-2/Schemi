/*
 * vanLeer2Limiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vanLeer2Limiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"
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
		const vector & gradient) const noexcept
{
	vector vanLeer2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), vanLeer2.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->vanLeer2LimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(vanLeer2()) * std::get<0>(gradient()), std::get<
			1>(vanLeer2()) * std::get<1>(gradient()), std::get<2>(vanLeer2())
			* std::get<2>(gradient()) };
}

schemi::tensor schemi::vanLeer2Limiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor vanLeer2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), vanLeer2.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->vanLeer2LimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(vanLeer2()) * std::get<0>(gradient()), std::get<
			1>(vanLeer2()) * std::get<1>(gradient()), std::get<2>(vanLeer2())
			* std::get<2>(gradient()), std::get<3>(vanLeer2())
			* std::get<3>(gradient()), std::get<4>(vanLeer2())
			* std::get<4>(gradient()), std::get<5>(vanLeer2())
			* std::get<5>(gradient()), std::get<6>(vanLeer2())
			* std::get<6>(gradient()), std::get<7>(vanLeer2())
			* std::get<7>(gradient()), std::get<8>(vanLeer2())
			* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::vanLeer2Limiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
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

	std::transform(r().begin(), r().end(), xiR().begin(), vanLeer2.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->vanLeer2LimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(vanLeer2()) * std::get<0>(gradient()),
			std::get<1>(vanLeer2()) * std::get<1>(gradient()), std::get<2>(
					vanLeer2()) * std::get<2>(gradient()), std::get<3>(
					vanLeer2()) * std::get<3>(gradient()), std::get<4>(
					vanLeer2()) * std::get<4>(gradient()), std::get<5>(
					vanLeer2()) * std::get<5>(gradient()), std::get<6>(
					vanLeer2()) * std::get<6>(gradient()), std::get<7>(
					vanLeer2()) * std::get<7>(gradient()), std::get<8>(
					vanLeer2()) * std::get<8>(gradient()), std::get<9>(
					vanLeer2()) * std::get<9>(gradient()), std::get<10>(
					vanLeer2()) * std::get<10>(gradient()), std::get<11>(
					vanLeer2()) * std::get<11>(gradient()), std::get<12>(
					vanLeer2()) * std::get<12>(gradient()), std::get<13>(
					vanLeer2()) * std::get<13>(gradient()), std::get<14>(
					vanLeer2()) * std::get<14>(gradient()), std::get<15>(
					vanLeer2()) * std::get<15>(gradient()), std::get<16>(
					vanLeer2()) * std::get<16>(gradient()), std::get<17>(
					vanLeer2()) * std::get<17>(gradient()), std::get<18>(
					vanLeer2()) * std::get<18>(gradient()), std::get<19>(
					vanLeer2()) * std::get<19>(gradient()), std::get<20>(
					vanLeer2()) * std::get<20>(gradient()), std::get<21>(
					vanLeer2()) * std::get<21>(gradient()), std::get<22>(
					vanLeer2()) * std::get<22>(gradient()), std::get<23>(
					vanLeer2()) * std::get<23>(gradient()), std::get<24>(
					vanLeer2()) * std::get<24>(gradient()), std::get<25>(
					vanLeer2()) * std::get<25>(gradient()), std::get<26>(
					vanLeer2()) * std::get<26>(gradient()) };
}

schemi::vector schemi::vanLeer2Limiter::calculate3OLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector vanLeer2;

	std::transform(r().begin(), r().end(), vanLeer2.r().begin(),
			[this](const auto r_j) 
			{	return this->vanLeer2LimiterCalculation(r_j);});

	const auto beta = elementsDivision((vector(1) + 2 * r), vector(3));

	std::transform(vanLeer2().begin(), vanLeer2().end(), beta().begin(),
			vanLeer2.r().begin(), [](const auto limiter_j, const auto beta_j) 
			{	return std::max(std::min(limiter_j, beta_j),0.0);});

	return vector { std::get<0>(vanLeer2()) * std::get<0>(gradient()), std::get<
			1>(vanLeer2()) * std::get<1>(gradient()), std::get<2>(vanLeer2())
			* std::get<2>(gradient()) };
}

schemi::tensor schemi::vanLeer2Limiter::calculate3OLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor vanLeer2;

	std::transform(r().begin(), r().end(), vanLeer2.r().begin(),
			[this](const auto r_j) 
			{	return this->vanLeer2LimiterCalculation(r_j);});

	const auto beta = elementsDivision((tensor(1) + 2 * r), tensor(3));

	std::transform(vanLeer2().begin(), vanLeer2().end(), beta().begin(),
			vanLeer2.r().begin(), [](const auto limiter_j, const auto beta_j) 
			{	return std::max(std::min(limiter_j, beta_j),0.0);});

	return tensor { std::get<0>(vanLeer2()) * std::get<0>(gradient()), std::get<
			1>(vanLeer2()) * std::get<1>(gradient()), std::get<2>(vanLeer2())
			* std::get<2>(gradient()), std::get<3>(vanLeer2())
			* std::get<3>(gradient()), std::get<4>(vanLeer2())
			* std::get<4>(gradient()), std::get<5>(vanLeer2())
			* std::get<5>(gradient()), std::get<6>(vanLeer2())
			* std::get<6>(gradient()), std::get<7>(vanLeer2())
			* std::get<7>(gradient()), std::get<8>(vanLeer2())
			* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::vanLeer2Limiter::calculate3OLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 vanLeer2;

	std::transform(r().begin(), r().end(), vanLeer2.r().begin(),
			[this](const auto r_j) 
			{	return this->vanLeer2LimiterCalculation(r_j);});

	const auto beta = elementsDivision((tensor3(1) + 2 * r), tensor3(3));

	std::transform(vanLeer2().begin(), vanLeer2().end(), beta().begin(),
			vanLeer2.r().begin(), [](const auto limiter_j, const auto beta_j) 
			{	return std::max(std::min(limiter_j, beta_j),0.0);});

	return tensor3 { std::get<0>(vanLeer2()) * std::get<0>(gradient()),
			std::get<1>(vanLeer2()) * std::get<1>(gradient()), std::get<2>(
					vanLeer2()) * std::get<2>(gradient()), std::get<3>(
					vanLeer2()) * std::get<3>(gradient()), std::get<4>(
					vanLeer2()) * std::get<4>(gradient()), std::get<5>(
					vanLeer2()) * std::get<5>(gradient()), std::get<6>(
					vanLeer2()) * std::get<6>(gradient()), std::get<7>(
					vanLeer2()) * std::get<7>(gradient()), std::get<8>(
					vanLeer2()) * std::get<8>(gradient()), std::get<9>(
					vanLeer2()) * std::get<9>(gradient()), std::get<10>(
					vanLeer2()) * std::get<10>(gradient()), std::get<11>(
					vanLeer2()) * std::get<11>(gradient()), std::get<12>(
					vanLeer2()) * std::get<12>(gradient()), std::get<13>(
					vanLeer2()) * std::get<13>(gradient()), std::get<14>(
					vanLeer2()) * std::get<14>(gradient()), std::get<15>(
					vanLeer2()) * std::get<15>(gradient()), std::get<16>(
					vanLeer2()) * std::get<16>(gradient()), std::get<17>(
					vanLeer2()) * std::get<17>(gradient()), std::get<18>(
					vanLeer2()) * std::get<18>(gradient()), std::get<19>(
					vanLeer2()) * std::get<19>(gradient()), std::get<20>(
					vanLeer2()) * std::get<20>(gradient()), std::get<21>(
					vanLeer2()) * std::get<21>(gradient()), std::get<22>(
					vanLeer2()) * std::get<22>(gradient()), std::get<23>(
					vanLeer2()) * std::get<23>(gradient()), std::get<24>(
					vanLeer2()) * std::get<24>(gradient()), std::get<25>(
					vanLeer2()) * std::get<25>(gradient()), std::get<26>(
					vanLeer2()) * std::get<26>(gradient()) };
}
