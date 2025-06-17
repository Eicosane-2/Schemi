/*
 * vanAlbada2Limiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vanAlbada2Limiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"
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
		const vector & gradient) const noexcept
{
	vector vanAlbada2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(),
			vanAlbada2.r().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->vanAlbada2LimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(vanAlbada2()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::vanAlbada2Limiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor vanAlbada2, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(),
			vanAlbada2.r().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->vanAlbada2LimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(vanAlbada2()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradient()), std::get<3>(
					vanAlbada2()) * std::get<3>(gradient()), std::get<4>(
					vanAlbada2()) * std::get<4>(gradient()), std::get<5>(
					vanAlbada2()) * std::get<5>(gradient()), std::get<6>(
					vanAlbada2()) * std::get<6>(gradient()), std::get<7>(
					vanAlbada2()) * std::get<7>(gradient()), std::get<8>(
					vanAlbada2()) * std::get<8>(gradient()) };
}

schemi::tensor3 schemi::vanAlbada2Limiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
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

	std::transform(r().begin(), r().end(), xiR().begin(),
			vanAlbada2.r().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->vanAlbada2LimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(vanAlbada2()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradient()), std::get<3>(
					vanAlbada2()) * std::get<3>(gradient()), std::get<4>(
					vanAlbada2()) * std::get<4>(gradient()), std::get<5>(
					vanAlbada2()) * std::get<5>(gradient()), std::get<6>(
					vanAlbada2()) * std::get<6>(gradient()), std::get<7>(
					vanAlbada2()) * std::get<7>(gradient()), std::get<8>(
					vanAlbada2()) * std::get<8>(gradient()), std::get<9>(
					vanAlbada2()) * std::get<9>(gradient()), std::get<10>(
					vanAlbada2()) * std::get<10>(gradient()), std::get<11>(
					vanAlbada2()) * std::get<11>(gradient()), std::get<12>(
					vanAlbada2()) * std::get<12>(gradient()), std::get<13>(
					vanAlbada2()) * std::get<13>(gradient()), std::get<14>(
					vanAlbada2()) * std::get<14>(gradient()), std::get<15>(
					vanAlbada2()) * std::get<15>(gradient()), std::get<16>(
					vanAlbada2()) * std::get<16>(gradient()), std::get<17>(
					vanAlbada2()) * std::get<17>(gradient()), std::get<18>(
					vanAlbada2()) * std::get<18>(gradient()), std::get<19>(
					vanAlbada2()) * std::get<19>(gradient()), std::get<20>(
					vanAlbada2()) * std::get<20>(gradient()), std::get<21>(
					vanAlbada2()) * std::get<21>(gradient()), std::get<22>(
					vanAlbada2()) * std::get<22>(gradient()), std::get<23>(
					vanAlbada2()) * std::get<23>(gradient()), std::get<24>(
					vanAlbada2()) * std::get<24>(gradient()), std::get<25>(
					vanAlbada2()) * std::get<25>(gradient()), std::get<26>(
					vanAlbada2()) * std::get<26>(gradient()) };
}

schemi::vector schemi::vanAlbada2Limiter::calculateNoRSLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector vanAlbada2;

	std::transform(r().begin(), r().end(), vanAlbada2.r().begin(),
			[this](const auto r_j) 
			{	return this->vanAlbada2LimiterCalculation(r_j);});

	return vector { std::get<0>(vanAlbada2()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::vanAlbada2Limiter::calculateNoRSLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor vanAlbada2;

	std::transform(r().begin(), r().end(), vanAlbada2.r().begin(),
			[this](const auto r_j) 
			{	return this->vanAlbada2LimiterCalculation(r_j);});

	return tensor { std::get<0>(vanAlbada2()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradient()), std::get<3>(
					vanAlbada2()) * std::get<3>(gradient()), std::get<4>(
					vanAlbada2()) * std::get<4>(gradient()), std::get<5>(
					vanAlbada2()) * std::get<5>(gradient()), std::get<6>(
					vanAlbada2()) * std::get<6>(gradient()), std::get<7>(
					vanAlbada2()) * std::get<7>(gradient()), std::get<8>(
					vanAlbada2()) * std::get<8>(gradient()) };
}

schemi::tensor3 schemi::vanAlbada2Limiter::calculateNoRSLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 vanAlbada2;

	std::transform(r().begin(), r().end(), vanAlbada2.r().begin(),
			[this](const auto r_j) 
			{	return this->vanAlbada2LimiterCalculation(r_j);});

	return tensor3 { std::get<0>(vanAlbada2()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada2()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada2()) * std::get<2>(gradient()), std::get<3>(
					vanAlbada2()) * std::get<3>(gradient()), std::get<4>(
					vanAlbada2()) * std::get<4>(gradient()), std::get<5>(
					vanAlbada2()) * std::get<5>(gradient()), std::get<6>(
					vanAlbada2()) * std::get<6>(gradient()), std::get<7>(
					vanAlbada2()) * std::get<7>(gradient()), std::get<8>(
					vanAlbada2()) * std::get<8>(gradient()), std::get<9>(
					vanAlbada2()) * std::get<9>(gradient()), std::get<10>(
					vanAlbada2()) * std::get<10>(gradient()), std::get<11>(
					vanAlbada2()) * std::get<11>(gradient()), std::get<12>(
					vanAlbada2()) * std::get<12>(gradient()), std::get<13>(
					vanAlbada2()) * std::get<13>(gradient()), std::get<14>(
					vanAlbada2()) * std::get<14>(gradient()), std::get<15>(
					vanAlbada2()) * std::get<15>(gradient()), std::get<16>(
					vanAlbada2()) * std::get<16>(gradient()), std::get<17>(
					vanAlbada2()) * std::get<17>(gradient()), std::get<18>(
					vanAlbada2()) * std::get<18>(gradient()), std::get<19>(
					vanAlbada2()) * std::get<19>(gradient()), std::get<20>(
					vanAlbada2()) * std::get<20>(gradient()), std::get<21>(
					vanAlbada2()) * std::get<21>(gradient()), std::get<22>(
					vanAlbada2()) * std::get<22>(gradient()), std::get<23>(
					vanAlbada2()) * std::get<23>(gradient()), std::get<24>(
					vanAlbada2()) * std::get<24>(gradient()), std::get<25>(
					vanAlbada2()) * std::get<25>(gradient()), std::get<26>(
					vanAlbada2()) * std::get<26>(gradient()) };
}
