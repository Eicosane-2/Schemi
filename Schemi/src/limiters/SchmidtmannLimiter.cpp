/*
 * SchmidtmannLimiter.cpp
 *
 *  Created on: 2025/05/24
 *      Author: Maxim Boldyrev
 */

#include "SchmidtmannLimiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"
#include "intExpPow.hpp"

schemi::scalar schemi::SchmidtmannLimiter::SchmidtmannLimiterCalculation(
		const scalar r, const scalar xiR) const noexcept
{
	scalar Schmidtmann;

	const auto phiO3 = (2 + r) / 3;

	Schmidtmann = std::min(
			std::max(static_cast<scalar>(0),
					std::min(phiO3, std::max(-alpha * r, std::min( { beta * r,
							phiO3, gamma })))), 2 * xiR);

	return Schmidtmann;
}

schemi::scalar schemi::SchmidtmannLimiter::SchmidtmannLimiterCalculation(
		const scalar r) const noexcept
{
	scalar Schmidtmann;

	const auto phiO3 = (2 + r) / 3;

	Schmidtmann = std::max(static_cast<scalar>(0),
			std::min(phiO3, std::max(-alpha * r, std::min( { beta * r, phiO3,
					gamma }))));

	return Schmidtmann;
}

schemi::vector schemi::SchmidtmannLimiter::calculate(const vector & r,
		const vector & gradient) const noexcept
{
	vector Schmidtmann, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			Schmidtmann.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->SchmidtmannLimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(Schmidtmann()) * std::get<0>(gradient()),
			std::get<1>(Schmidtmann()) * std::get<1>(gradient()), std::get<2>(
					Schmidtmann()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::SchmidtmannLimiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor Schmidtmann, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			Schmidtmann.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->SchmidtmannLimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(Schmidtmann()) * std::get<0>(gradient()),
			std::get<1>(Schmidtmann()) * std::get<1>(gradient()), std::get<2>(
					Schmidtmann()) * std::get<2>(gradient()), std::get<3>(
					Schmidtmann()) * std::get<3>(gradient()), std::get<4>(
					Schmidtmann()) * std::get<4>(gradient()), std::get<5>(
					Schmidtmann()) * std::get<5>(gradient()), std::get<6>(
					Schmidtmann()) * std::get<6>(gradient()), std::get<7>(
					Schmidtmann()) * std::get<7>(gradient()), std::get<8>(
					Schmidtmann()) * std::get<8>(gradient()) };
}

schemi::tensor3 schemi::SchmidtmannLimiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 Schmidtmann, xiR { 2 / (1 + std::get<0>(r())), 2
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
			Schmidtmann.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->SchmidtmannLimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(Schmidtmann()) * std::get<0>(gradient()),
			std::get<1>(Schmidtmann()) * std::get<1>(gradient()), std::get<2>(
					Schmidtmann()) * std::get<2>(gradient()), std::get<3>(
					Schmidtmann()) * std::get<3>(gradient()), std::get<4>(
					Schmidtmann()) * std::get<4>(gradient()), std::get<5>(
					Schmidtmann()) * std::get<5>(gradient()), std::get<6>(
					Schmidtmann()) * std::get<6>(gradient()), std::get<7>(
					Schmidtmann()) * std::get<7>(gradient()), std::get<8>(
					Schmidtmann()) * std::get<8>(gradient()), std::get<9>(
					Schmidtmann()) * std::get<9>(gradient()), std::get<10>(
					Schmidtmann()) * std::get<10>(gradient()), std::get<11>(
					Schmidtmann()) * std::get<11>(gradient()), std::get<12>(
					Schmidtmann()) * std::get<12>(gradient()), std::get<13>(
					Schmidtmann()) * std::get<13>(gradient()), std::get<14>(
					Schmidtmann()) * std::get<14>(gradient()), std::get<15>(
					Schmidtmann()) * std::get<15>(gradient()), std::get<16>(
					Schmidtmann()) * std::get<16>(gradient()), std::get<17>(
					Schmidtmann()) * std::get<17>(gradient()), std::get<18>(
					Schmidtmann()) * std::get<18>(gradient()), std::get<19>(
					Schmidtmann()) * std::get<19>(gradient()), std::get<20>(
					Schmidtmann()) * std::get<20>(gradient()), std::get<21>(
					Schmidtmann()) * std::get<21>(gradient()), std::get<22>(
					Schmidtmann()) * std::get<22>(gradient()), std::get<23>(
					Schmidtmann()) * std::get<23>(gradient()), std::get<24>(
					Schmidtmann()) * std::get<24>(gradient()), std::get<25>(
					Schmidtmann()) * std::get<25>(gradient()), std::get<26>(
					Schmidtmann()) * std::get<26>(gradient()) };
}

schemi::vector schemi::SchmidtmannLimiter::calculateNoRSLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector Schmidtmann;

	std::transform(r().cbegin(), r().cend(), Schmidtmann.wr().begin(),
			[this](const auto r_j) 
			{	return this->SchmidtmannLimiterCalculation(r_j);});

	return vector { std::get<0>(Schmidtmann()) * std::get<0>(gradient()),
			std::get<1>(Schmidtmann()) * std::get<1>(gradient()), std::get<2>(
					Schmidtmann()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::SchmidtmannLimiter::calculateNoRSLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor Schmidtmann;

	std::transform(r().cbegin(), r().cend(), Schmidtmann.wr().begin(),
			[this](const auto r_j) 
			{	return this->SchmidtmannLimiterCalculation(r_j);});

	return tensor { std::get<0>(Schmidtmann()) * std::get<0>(gradient()),
			std::get<1>(Schmidtmann()) * std::get<1>(gradient()), std::get<2>(
					Schmidtmann()) * std::get<2>(gradient()), std::get<3>(
					Schmidtmann()) * std::get<3>(gradient()), std::get<4>(
					Schmidtmann()) * std::get<4>(gradient()), std::get<5>(
					Schmidtmann()) * std::get<5>(gradient()), std::get<6>(
					Schmidtmann()) * std::get<6>(gradient()), std::get<7>(
					Schmidtmann()) * std::get<7>(gradient()), std::get<8>(
					Schmidtmann()) * std::get<8>(gradient()) };
}

schemi::tensor3 schemi::SchmidtmannLimiter::calculateNoRSLimit(
		const tensor3 & r, const tensor3 & gradient) const noexcept
{
	tensor3 Schmidtmann;

	std::transform(r().cbegin(), r().cend(), Schmidtmann.wr().begin(),
			[this](const auto r_j) 
			{	return this->SchmidtmannLimiterCalculation(r_j);});

	return tensor3 { std::get<0>(Schmidtmann()) * std::get<0>(gradient()),
			std::get<1>(Schmidtmann()) * std::get<1>(gradient()), std::get<2>(
					Schmidtmann()) * std::get<2>(gradient()), std::get<3>(
					Schmidtmann()) * std::get<3>(gradient()), std::get<4>(
					Schmidtmann()) * std::get<4>(gradient()), std::get<5>(
					Schmidtmann()) * std::get<5>(gradient()), std::get<6>(
					Schmidtmann()) * std::get<6>(gradient()), std::get<7>(
					Schmidtmann()) * std::get<7>(gradient()), std::get<8>(
					Schmidtmann()) * std::get<8>(gradient()), std::get<9>(
					Schmidtmann()) * std::get<9>(gradient()), std::get<10>(
					Schmidtmann()) * std::get<10>(gradient()), std::get<11>(
					Schmidtmann()) * std::get<11>(gradient()), std::get<12>(
					Schmidtmann()) * std::get<12>(gradient()), std::get<13>(
					Schmidtmann()) * std::get<13>(gradient()), std::get<14>(
					Schmidtmann()) * std::get<14>(gradient()), std::get<15>(
					Schmidtmann()) * std::get<15>(gradient()), std::get<16>(
					Schmidtmann()) * std::get<16>(gradient()), std::get<17>(
					Schmidtmann()) * std::get<17>(gradient()), std::get<18>(
					Schmidtmann()) * std::get<18>(gradient()), std::get<19>(
					Schmidtmann()) * std::get<19>(gradient()), std::get<20>(
					Schmidtmann()) * std::get<20>(gradient()), std::get<21>(
					Schmidtmann()) * std::get<21>(gradient()), std::get<22>(
					Schmidtmann()) * std::get<22>(gradient()), std::get<23>(
					Schmidtmann()) * std::get<23>(gradient()), std::get<24>(
					Schmidtmann()) * std::get<24>(gradient()), std::get<25>(
					Schmidtmann()) * std::get<25>(gradient()), std::get<26>(
					Schmidtmann()) * std::get<26>(gradient()) };
}
