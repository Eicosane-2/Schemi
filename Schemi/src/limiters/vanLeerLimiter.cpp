/*
 * vanLeerLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vanLeerLimiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"

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
		const vector & gradient) const noexcept
{
	vector vanLeer, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			vanLeer.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->vanLeerLimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(vanLeer()) * std::get<0>(gradient()),
			std::get<1>(vanLeer()) * std::get<1>(gradient()), std::get<2>(
					vanLeer()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::vanLeerLimiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor vanLeer, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			vanLeer.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->vanLeerLimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(vanLeer()) * std::get<0>(gradient()),
			std::get<1>(vanLeer()) * std::get<1>(gradient()), std::get<2>(
					vanLeer()) * std::get<2>(gradient()), std::get<3>(vanLeer())
					* std::get<3>(gradient()), std::get<4>(vanLeer())
					* std::get<4>(gradient()), std::get<5>(vanLeer())
					* std::get<5>(gradient()), std::get<6>(vanLeer())
					* std::get<6>(gradient()), std::get<7>(vanLeer())
					* std::get<7>(gradient()), std::get<8>(vanLeer())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::vanLeerLimiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
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

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			vanLeer.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->vanLeerLimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(vanLeer()) * std::get<0>(gradient()), std::get<
			1>(vanLeer()) * std::get<1>(gradient()), std::get<2>(vanLeer())
			* std::get<2>(gradient()), std::get<3>(vanLeer())
			* std::get<3>(gradient()), std::get<4>(vanLeer())
			* std::get<4>(gradient()), std::get<5>(vanLeer())
			* std::get<5>(gradient()), std::get<6>(vanLeer())
			* std::get<6>(gradient()), std::get<7>(vanLeer())
			* std::get<7>(gradient()), std::get<8>(vanLeer())
			* std::get<8>(gradient()), std::get<9>(vanLeer())
			* std::get<9>(gradient()), std::get<10>(vanLeer())
			* std::get<10>(gradient()), std::get<11>(vanLeer())
			* std::get<11>(gradient()), std::get<12>(vanLeer())
			* std::get<12>(gradient()), std::get<13>(vanLeer())
			* std::get<13>(gradient()), std::get<14>(vanLeer())
			* std::get<14>(gradient()), std::get<15>(vanLeer())
			* std::get<15>(gradient()), std::get<16>(vanLeer())
			* std::get<16>(gradient()), std::get<17>(vanLeer())
			* std::get<17>(gradient()), std::get<18>(vanLeer())
			* std::get<18>(gradient()), std::get<19>(vanLeer())
			* std::get<19>(gradient()), std::get<20>(vanLeer())
			* std::get<20>(gradient()), std::get<21>(vanLeer())
			* std::get<21>(gradient()), std::get<22>(vanLeer())
			* std::get<22>(gradient()), std::get<23>(vanLeer())
			* std::get<23>(gradient()), std::get<24>(vanLeer())
			* std::get<24>(gradient()), std::get<25>(vanLeer())
			* std::get<25>(gradient()), std::get<26>(vanLeer())
			* std::get<26>(gradient()) };
}

schemi::vector schemi::vanLeerLimiter::calculateNoRSLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector vanLeer;

	std::transform(r().cbegin(), r().cend(), vanLeer.wr().begin(),
			[this](const auto r_j) 
			{	return this->vanLeerLimiterCalculation(r_j);});

	return vector { std::get<0>(vanLeer()) * std::get<0>(gradient()),
			std::get<1>(vanLeer()) * std::get<1>(gradient()), std::get<2>(
					vanLeer()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::vanLeerLimiter::calculateNoRSLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor vanLeer;

	std::transform(r().cbegin(), r().cend(), vanLeer.wr().begin(),
			[this](const auto r_j) 
			{	return this->vanLeerLimiterCalculation(r_j);});

	return tensor { std::get<0>(vanLeer()) * std::get<0>(gradient()),
			std::get<1>(vanLeer()) * std::get<1>(gradient()), std::get<2>(
					vanLeer()) * std::get<2>(gradient()), std::get<3>(vanLeer())
					* std::get<3>(gradient()), std::get<4>(vanLeer())
					* std::get<4>(gradient()), std::get<5>(vanLeer())
					* std::get<5>(gradient()), std::get<6>(vanLeer())
					* std::get<6>(gradient()), std::get<7>(vanLeer())
					* std::get<7>(gradient()), std::get<8>(vanLeer())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::vanLeerLimiter::calculateNoRSLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 vanLeer;

	std::transform(r().cbegin(), r().cend(), vanLeer.wr().begin(),
			[this](const auto r_j) 
			{	return this->vanLeerLimiterCalculation(r_j);});

	return tensor3 { std::get<0>(vanLeer()) * std::get<0>(gradient()), std::get<
			1>(vanLeer()) * std::get<1>(gradient()), std::get<2>(vanLeer())
			* std::get<2>(gradient()), std::get<3>(vanLeer())
			* std::get<3>(gradient()), std::get<4>(vanLeer())
			* std::get<4>(gradient()), std::get<5>(vanLeer())
			* std::get<5>(gradient()), std::get<6>(vanLeer())
			* std::get<6>(gradient()), std::get<7>(vanLeer())
			* std::get<7>(gradient()), std::get<8>(vanLeer())
			* std::get<8>(gradient()), std::get<9>(vanLeer())
			* std::get<9>(gradient()), std::get<10>(vanLeer())
			* std::get<10>(gradient()), std::get<11>(vanLeer())
			* std::get<11>(gradient()), std::get<12>(vanLeer())
			* std::get<12>(gradient()), std::get<13>(vanLeer())
			* std::get<13>(gradient()), std::get<14>(vanLeer())
			* std::get<14>(gradient()), std::get<15>(vanLeer())
			* std::get<15>(gradient()), std::get<16>(vanLeer())
			* std::get<16>(gradient()), std::get<17>(vanLeer())
			* std::get<17>(gradient()), std::get<18>(vanLeer())
			* std::get<18>(gradient()), std::get<19>(vanLeer())
			* std::get<19>(gradient()), std::get<20>(vanLeer())
			* std::get<20>(gradient()), std::get<21>(vanLeer())
			* std::get<21>(gradient()), std::get<22>(vanLeer())
			* std::get<22>(gradient()), std::get<23>(vanLeer())
			* std::get<23>(gradient()), std::get<24>(vanLeer())
			* std::get<24>(gradient()), std::get<25>(vanLeer())
			* std::get<25>(gradient()), std::get<26>(vanLeer())
			* std::get<26>(gradient()) };
}
