/*
 * superbeeLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "superbeeLimiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"

schemi::scalar schemi::superbeeLimiter::superbeeLimiterCalculation(
		const scalar r, const scalar xiR) const noexcept
{
	scalar superbee;

	if (r <= 0.)
		superbee = 0;
	else if (r <= 0.5)
		superbee = 2 * r;
	else if (r <= 1.)
		superbee = r;
	else
		superbee = std::min( { r, static_cast<scalar>(2), xiR });

	return superbee;
}

schemi::scalar schemi::superbeeLimiter::superbeeLimiterCalculation(
		const scalar r) const noexcept
{
	scalar superbee;

	if (r <= 0.)
		superbee = 0;
	else if (r <= 0.5)
		superbee = 2 * r;
	else if (r <= 1.)
		superbee = r;
	else
		superbee = std::min(r, static_cast<scalar>(2));

	return superbee;
}

schemi::vector schemi::superbeeLimiter::calculate(const vector & r,
		const vector & gradient) const noexcept
{
	vector superbee, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), superbee.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->superbeeLimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(superbee()) * std::get<0>(gradient()), std::get<
			1>(superbee()) * std::get<1>(gradient()), std::get<2>(superbee())
			* std::get<2>(gradient()) };
}

schemi::tensor schemi::superbeeLimiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor superbee, xiR { 2 / (1 + std::get<0>(r())), 2
			/ (1 + std::get<1>(r())), 2 / (1 + std::get<2>(r())), 2
			/ (1 + std::get<3>(r())), 2 / (1 + std::get<4>(r())), 2
			/ (1 + std::get<5>(r())), 2 / (1 + std::get<6>(r())), 2
			/ (1 + std::get<7>(r())), 2 / (1 + std::get<8>(r())) };

	std::transform(r().begin(), r().end(), xiR().begin(), superbee.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->superbeeLimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(superbee()) * std::get<0>(gradient()), std::get<
			1>(superbee()) * std::get<1>(gradient()), std::get<2>(superbee())
			* std::get<2>(gradient()), std::get<3>(superbee())
			* std::get<3>(gradient()), std::get<4>(superbee())
			* std::get<4>(gradient()), std::get<5>(superbee())
			* std::get<5>(gradient()), std::get<6>(superbee())
			* std::get<6>(gradient()), std::get<7>(superbee())
			* std::get<7>(gradient()), std::get<8>(superbee())
			* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::superbeeLimiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 superbee, xiR { 2 / (1 + std::get<0>(r())), 2
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

	std::transform(r().begin(), r().end(), xiR().begin(), superbee.r().begin(),
			[this](const auto r_j, const auto xiR_j) 
			{	return this->superbeeLimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(superbee()) * std::get<0>(gradient()),
			std::get<1>(superbee()) * std::get<1>(gradient()), std::get<2>(
					superbee()) * std::get<2>(gradient()), std::get<3>(
					superbee()) * std::get<3>(gradient()), std::get<4>(
					superbee()) * std::get<4>(gradient()), std::get<5>(
					superbee()) * std::get<5>(gradient()), std::get<6>(
					superbee()) * std::get<6>(gradient()), std::get<7>(
					superbee()) * std::get<7>(gradient()), std::get<8>(
					superbee()) * std::get<8>(gradient()), std::get<9>(
					superbee()) * std::get<9>(gradient()), std::get<10>(
					superbee()) * std::get<10>(gradient()), std::get<11>(
					superbee()) * std::get<11>(gradient()), std::get<12>(
					superbee()) * std::get<12>(gradient()), std::get<13>(
					superbee()) * std::get<13>(gradient()), std::get<14>(
					superbee()) * std::get<14>(gradient()), std::get<15>(
					superbee()) * std::get<15>(gradient()), std::get<16>(
					superbee()) * std::get<16>(gradient()), std::get<17>(
					superbee()) * std::get<17>(gradient()), std::get<18>(
					superbee()) * std::get<18>(gradient()), std::get<19>(
					superbee()) * std::get<19>(gradient()), std::get<20>(
					superbee()) * std::get<20>(gradient()), std::get<21>(
					superbee()) * std::get<21>(gradient()), std::get<22>(
					superbee()) * std::get<22>(gradient()), std::get<23>(
					superbee()) * std::get<23>(gradient()), std::get<24>(
					superbee()) * std::get<24>(gradient()), std::get<25>(
					superbee()) * std::get<25>(gradient()), std::get<26>(
					superbee()) * std::get<26>(gradient()) };
}

schemi::vector schemi::superbeeLimiter::calculateNoRSLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector superbee;

	std::transform(r().begin(), r().end(), superbee.r().begin(),
			[this](const auto r_j) 
			{	return this->superbeeLimiterCalculation(r_j);});

	return vector { std::get<0>(superbee()) * std::get<0>(gradient()), std::get<
			1>(superbee()) * std::get<1>(gradient()), std::get<2>(superbee())
			* std::get<2>(gradient()) };
}

schemi::tensor schemi::superbeeLimiter::calculateNoRSLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor superbee;

	std::transform(r().begin(), r().end(), superbee.r().begin(),
			[this](const auto r_j) 
			{	return this->superbeeLimiterCalculation(r_j);});

	return tensor { std::get<0>(superbee()) * std::get<0>(gradient()), std::get<
			1>(superbee()) * std::get<1>(gradient()), std::get<2>(superbee())
			* std::get<2>(gradient()), std::get<3>(superbee())
			* std::get<3>(gradient()), std::get<4>(superbee())
			* std::get<4>(gradient()), std::get<5>(superbee())
			* std::get<5>(gradient()), std::get<6>(superbee())
			* std::get<6>(gradient()), std::get<7>(superbee())
			* std::get<7>(gradient()), std::get<8>(superbee())
			* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::superbeeLimiter::calculateNoRSLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 superbee;

	std::transform(r().begin(), r().end(), superbee.r().begin(),
			[this](const auto r_j) 
			{	return this->superbeeLimiterCalculation(r_j);});

	return tensor3 { std::get<0>(superbee()) * std::get<0>(gradient()),
			std::get<1>(superbee()) * std::get<1>(gradient()), std::get<2>(
					superbee()) * std::get<2>(gradient()), std::get<3>(
					superbee()) * std::get<3>(gradient()), std::get<4>(
					superbee()) * std::get<4>(gradient()), std::get<5>(
					superbee()) * std::get<5>(gradient()), std::get<6>(
					superbee()) * std::get<6>(gradient()), std::get<7>(
					superbee()) * std::get<7>(gradient()), std::get<8>(
					superbee()) * std::get<8>(gradient()), std::get<9>(
					superbee()) * std::get<9>(gradient()), std::get<10>(
					superbee()) * std::get<10>(gradient()), std::get<11>(
					superbee()) * std::get<11>(gradient()), std::get<12>(
					superbee()) * std::get<12>(gradient()), std::get<13>(
					superbee()) * std::get<13>(gradient()), std::get<14>(
					superbee()) * std::get<14>(gradient()), std::get<15>(
					superbee()) * std::get<15>(gradient()), std::get<16>(
					superbee()) * std::get<16>(gradient()), std::get<17>(
					superbee()) * std::get<17>(gradient()), std::get<18>(
					superbee()) * std::get<18>(gradient()), std::get<19>(
					superbee()) * std::get<19>(gradient()), std::get<20>(
					superbee()) * std::get<20>(gradient()), std::get<21>(
					superbee()) * std::get<21>(gradient()), std::get<22>(
					superbee()) * std::get<22>(gradient()), std::get<23>(
					superbee()) * std::get<23>(gradient()), std::get<24>(
					superbee()) * std::get<24>(gradient()), std::get<25>(
					superbee()) * std::get<25>(gradient()), std::get<26>(
					superbee()) * std::get<26>(gradient()) };
}
