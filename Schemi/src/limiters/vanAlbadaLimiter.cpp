/*
 * vanAlbadaLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vanAlbadaLimiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"
#include "intExpPow.hpp"

schemi::scalar schemi::vanAlbadaLimiter::vanAlbadaLimiterCalculation(
		const scalar r, const scalar xiR) const noexcept
{
	scalar vanAlbada;

	if (r <= 0.)
		vanAlbada = 0;
	else
		vanAlbada = std::min((pow<scalar, 2>(r) + r) / (1 + pow<scalar, 2>(r)),
				xiR);

	return vanAlbada;
}

schemi::scalar schemi::vanAlbadaLimiter::vanAlbadaLimiterCalculation(
		const scalar r) const noexcept
{
	scalar vanAlbada;

	if (r <= 0.)
		vanAlbada = 0;
	else
		vanAlbada = (pow<scalar, 2>(r) + r) / (1 + pow<scalar, 2>(r));

	return vanAlbada;
}

schemi::vector schemi::vanAlbadaLimiter::calculate(const vector & r,
		const vector & gradient) const noexcept
{
	vector vanAlbada, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), vanAlbada.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->vanAlbadaLimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(vanAlbada()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::vanAlbadaLimiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor vanAlbada, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), vanAlbada.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->vanAlbadaLimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(vanAlbada()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada()) * std::get<2>(gradient()), std::get<3>(
					vanAlbada()) * std::get<3>(gradient()), std::get<4>(
					vanAlbada()) * std::get<4>(gradient()), std::get<5>(
					vanAlbada()) * std::get<5>(gradient()), std::get<6>(
					vanAlbada()) * std::get<6>(gradient()), std::get<7>(
					vanAlbada()) * std::get<7>(gradient()), std::get<8>(
					vanAlbada()) * std::get<8>(gradient()) };
}

schemi::tensor3 schemi::vanAlbadaLimiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 vanAlbada, xiR { 2 / (1 + std::get<0>(r())), 2
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

	std::transform(r().begin(), r().end(), xiR().begin(), vanAlbada.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->vanAlbadaLimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(vanAlbada()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada()) * std::get<2>(gradient()), std::get<3>(
					vanAlbada()) * std::get<3>(gradient()), std::get<4>(
					vanAlbada()) * std::get<4>(gradient()), std::get<5>(
					vanAlbada()) * std::get<5>(gradient()), std::get<6>(
					vanAlbada()) * std::get<6>(gradient()), std::get<7>(
					vanAlbada()) * std::get<7>(gradient()), std::get<8>(
					vanAlbada()) * std::get<8>(gradient()), std::get<9>(
					vanAlbada()) * std::get<9>(gradient()), std::get<10>(
					vanAlbada()) * std::get<10>(gradient()), std::get<11>(
					vanAlbada()) * std::get<11>(gradient()), std::get<12>(
					vanAlbada()) * std::get<12>(gradient()), std::get<13>(
					vanAlbada()) * std::get<13>(gradient()), std::get<14>(
					vanAlbada()) * std::get<14>(gradient()), std::get<15>(
					vanAlbada()) * std::get<15>(gradient()), std::get<16>(
					vanAlbada()) * std::get<16>(gradient()), std::get<17>(
					vanAlbada()) * std::get<17>(gradient()), std::get<18>(
					vanAlbada()) * std::get<18>(gradient()), std::get<19>(
					vanAlbada()) * std::get<19>(gradient()), std::get<20>(
					vanAlbada()) * std::get<20>(gradient()), std::get<21>(
					vanAlbada()) * std::get<21>(gradient()), std::get<22>(
					vanAlbada()) * std::get<22>(gradient()), std::get<23>(
					vanAlbada()) * std::get<23>(gradient()), std::get<24>(
					vanAlbada()) * std::get<24>(gradient()), std::get<25>(
					vanAlbada()) * std::get<25>(gradient()), std::get<26>(
					vanAlbada()) * std::get<26>(gradient()) };
}

schemi::vector schemi::vanAlbadaLimiter::calculateNoRSLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector vanAlbada;

	std::transform(r().begin(), r().end(), vanAlbada.r().begin(),
			[this](const auto r_j) 
			{	return this->vanAlbadaLimiterCalculation(r_j);});

	return vector { std::get<0>(vanAlbada()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::vanAlbadaLimiter::calculateNoRSLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor vanAlbada;

	std::transform(r().begin(), r().end(), vanAlbada.r().begin(),
			[this](const auto r_j) 
			{	return this->vanAlbadaLimiterCalculation(r_j);});

	return tensor { std::get<0>(vanAlbada()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada()) * std::get<2>(gradient()), std::get<3>(
					vanAlbada()) * std::get<3>(gradient()), std::get<4>(
					vanAlbada()) * std::get<4>(gradient()), std::get<5>(
					vanAlbada()) * std::get<5>(gradient()), std::get<6>(
					vanAlbada()) * std::get<6>(gradient()), std::get<7>(
					vanAlbada()) * std::get<7>(gradient()), std::get<8>(
					vanAlbada()) * std::get<8>(gradient()) };
}

schemi::tensor3 schemi::vanAlbadaLimiter::calculateNoRSLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 vanAlbada;

	std::transform(r().begin(), r().end(), vanAlbada.r().begin(),
			[this](const auto r_j) 
			{	return this->vanAlbadaLimiterCalculation(r_j);});

	return tensor3 { std::get<0>(vanAlbada()) * std::get<0>(gradient()),
			std::get<1>(vanAlbada()) * std::get<1>(gradient()), std::get<2>(
					vanAlbada()) * std::get<2>(gradient()), std::get<3>(
					vanAlbada()) * std::get<3>(gradient()), std::get<4>(
					vanAlbada()) * std::get<4>(gradient()), std::get<5>(
					vanAlbada()) * std::get<5>(gradient()), std::get<6>(
					vanAlbada()) * std::get<6>(gradient()), std::get<7>(
					vanAlbada()) * std::get<7>(gradient()), std::get<8>(
					vanAlbada()) * std::get<8>(gradient()), std::get<9>(
					vanAlbada()) * std::get<9>(gradient()), std::get<10>(
					vanAlbada()) * std::get<10>(gradient()), std::get<11>(
					vanAlbada()) * std::get<11>(gradient()), std::get<12>(
					vanAlbada()) * std::get<12>(gradient()), std::get<13>(
					vanAlbada()) * std::get<13>(gradient()), std::get<14>(
					vanAlbada()) * std::get<14>(gradient()), std::get<15>(
					vanAlbada()) * std::get<15>(gradient()), std::get<16>(
					vanAlbada()) * std::get<16>(gradient()), std::get<17>(
					vanAlbada()) * std::get<17>(gradient()), std::get<18>(
					vanAlbada()) * std::get<18>(gradient()), std::get<19>(
					vanAlbada()) * std::get<19>(gradient()), std::get<20>(
					vanAlbada()) * std::get<20>(gradient()), std::get<21>(
					vanAlbada()) * std::get<21>(gradient()), std::get<22>(
					vanAlbada()) * std::get<22>(gradient()), std::get<23>(
					vanAlbada()) * std::get<23>(gradient()), std::get<24>(
					vanAlbada()) * std::get<24>(gradient()), std::get<25>(
					vanAlbada()) * std::get<25>(gradient()), std::get<26>(
					vanAlbada()) * std::get<26>(gradient()) };
}
