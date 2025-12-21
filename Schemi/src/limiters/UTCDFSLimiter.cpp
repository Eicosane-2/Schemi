/*
 * UTCDFSLimiter.cpp
 *
 *  Created on: 2025/05/23
 *      Author: Maxim Boldyrev
 */

#include "UTCDFSLimiter.hpp"

#include <algorithm>

#include "elementsProduct.hpp"
#include "intExpPow.hpp"

schemi::scalar schemi::UTCDFSLimiter::UTCDFSLimiterCalculation(const scalar r,
		const scalar xiR) const noexcept
{
	scalar UTCDFS;

	if (r < -1)
		UTCDFS = r * (r + 1) / (pow<scalar, 2>(r) + 1);
	else if ((r >= -1) && (r < 0))
		UTCDFS = 0;
	else if ((r >= 0) && (r < 0.5))
		UTCDFS = pow<scalar, 3>(r) - 2 * pow<scalar, 2>(r) + 2 * r;
	else if ((r >= 0.5) && (r < 2))
		UTCDFS = 0.75 * r + 0.25;
	else
		UTCDFS = std::min(
				(2 * pow<scalar, 2>(r) - 2 * r - 9.0 / 4.0)
						/ (pow<scalar, 2>(r) - r - 1), 2 * xiR);

	return UTCDFS;
}

schemi::scalar schemi::UTCDFSLimiter::UTCDFSLimiterCalculation(
		const scalar r) const noexcept
{
	scalar UTCDFS;

	if (r < -1)
		UTCDFS = r * (r + 1) / (pow<scalar, 2>(r) + 1);
	else if ((r >= -1) && (r < 0))
		UTCDFS = 0;
	else if ((r >= 0) && (r < 0.5))
		UTCDFS = pow<scalar, 3>(r) - 2 * pow<scalar, 2>(r) + 2 * r;
	else if ((r >= 0.5) && (r < 2))
		UTCDFS = 0.75 * r + 0.25;
	else
		UTCDFS = (2 * pow<scalar, 2>(r) - 2 * r - 9.0 / 4.0)
				/ (pow<scalar, 2>(r) - r - 1);

	return UTCDFS;
}

schemi::vector schemi::UTCDFSLimiter::calculate(const vector & r,
		const vector & gradient) const noexcept
{
	vector UTCDFS, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())) };

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			UTCDFS.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->UTCDFSLimiterCalculation(r_j, xiR_j);});

	return vector { std::get<0>(UTCDFS()) * std::get<0>(gradient()),
			std::get<1>(UTCDFS()) * std::get<1>(gradient()), std::get<2>(
					UTCDFS()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::UTCDFSLimiter::calculate(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor UTCDFS, xiR { 2 / (1 + std::get<0>(r())), 2 / (1 + std::get<1>(r())),
			2 / (1 + std::get<2>(r())), 2 / (1 + std::get<3>(r())), 2
					/ (1 + std::get<4>(r())), 2 / (1 + std::get<5>(r())), 2
					/ (1 + std::get<6>(r())), 2 / (1 + std::get<7>(r())), 2
					/ (1 + std::get<8>(r())) };

	std::transform(r().cbegin(), r().cend(), xiR().cbegin(),
			UTCDFS.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->UTCDFSLimiterCalculation(r_j, xiR_j);});

	return tensor { std::get<0>(UTCDFS()) * std::get<0>(gradient()),
			std::get<1>(UTCDFS()) * std::get<1>(gradient()), std::get<2>(
					UTCDFS()) * std::get<2>(gradient()), std::get<3>(UTCDFS())
					* std::get<3>(gradient()), std::get<4>(UTCDFS())
					* std::get<4>(gradient()), std::get<5>(UTCDFS())
					* std::get<5>(gradient()), std::get<6>(UTCDFS())
					* std::get<6>(gradient()), std::get<7>(UTCDFS())
					* std::get<7>(gradient()), std::get<8>(UTCDFS())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::UTCDFSLimiter::calculate(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 UTCDFS, xiR { 2 / (1 + std::get<0>(r())), 2
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
			UTCDFS.wr().begin(), [this](const auto r_j, const auto xiR_j) 
			{	return this->UTCDFSLimiterCalculation(r_j, xiR_j);});

	return tensor3 { std::get<0>(UTCDFS()) * std::get<0>(gradient()),
			std::get<1>(UTCDFS()) * std::get<1>(gradient()), std::get<2>(
					UTCDFS()) * std::get<2>(gradient()), std::get<3>(UTCDFS())
					* std::get<3>(gradient()), std::get<4>(UTCDFS())
					* std::get<4>(gradient()), std::get<5>(UTCDFS())
					* std::get<5>(gradient()), std::get<6>(UTCDFS())
					* std::get<6>(gradient()), std::get<7>(UTCDFS())
					* std::get<7>(gradient()), std::get<8>(UTCDFS())
					* std::get<8>(gradient()), std::get<9>(UTCDFS())
					* std::get<9>(gradient()), std::get<10>(UTCDFS())
					* std::get<10>(gradient()), std::get<11>(UTCDFS())
					* std::get<11>(gradient()), std::get<12>(UTCDFS())
					* std::get<12>(gradient()), std::get<13>(UTCDFS())
					* std::get<13>(gradient()), std::get<14>(UTCDFS())
					* std::get<14>(gradient()), std::get<15>(UTCDFS())
					* std::get<15>(gradient()), std::get<16>(UTCDFS())
					* std::get<16>(gradient()), std::get<17>(UTCDFS())
					* std::get<17>(gradient()), std::get<18>(UTCDFS())
					* std::get<18>(gradient()), std::get<19>(UTCDFS())
					* std::get<19>(gradient()), std::get<20>(UTCDFS())
					* std::get<20>(gradient()), std::get<21>(UTCDFS())
					* std::get<21>(gradient()), std::get<22>(UTCDFS())
					* std::get<22>(gradient()), std::get<23>(UTCDFS())
					* std::get<23>(gradient()), std::get<24>(UTCDFS())
					* std::get<24>(gradient()), std::get<25>(UTCDFS())
					* std::get<25>(gradient()), std::get<26>(UTCDFS())
					* std::get<26>(gradient()) };
}

schemi::vector schemi::UTCDFSLimiter::calculateNoRSLimit(const vector & r,
		const vector & gradient) const noexcept
{
	vector UTCDFS;

	std::transform(r().cbegin(), r().cend(), UTCDFS.wr().begin(),
			[this](const auto r_j) 
			{	return this->UTCDFSLimiterCalculation(r_j);});

	return vector { std::get<0>(UTCDFS()) * std::get<0>(gradient()),
			std::get<1>(UTCDFS()) * std::get<1>(gradient()), std::get<2>(
					UTCDFS()) * std::get<2>(gradient()) };
}

schemi::tensor schemi::UTCDFSLimiter::calculateNoRSLimit(const tensor & r,
		const tensor & gradient) const noexcept
{
	tensor UTCDFS;

	std::transform(r().cbegin(), r().cend(), UTCDFS.wr().begin(),
			[this](const auto r_j) 
			{	return this->UTCDFSLimiterCalculation(r_j);});

	return tensor { std::get<0>(UTCDFS()) * std::get<0>(gradient()),
			std::get<1>(UTCDFS()) * std::get<1>(gradient()), std::get<2>(
					UTCDFS()) * std::get<2>(gradient()), std::get<3>(UTCDFS())
					* std::get<3>(gradient()), std::get<4>(UTCDFS())
					* std::get<4>(gradient()), std::get<5>(UTCDFS())
					* std::get<5>(gradient()), std::get<6>(UTCDFS())
					* std::get<6>(gradient()), std::get<7>(UTCDFS())
					* std::get<7>(gradient()), std::get<8>(UTCDFS())
					* std::get<8>(gradient()) };
}

schemi::tensor3 schemi::UTCDFSLimiter::calculateNoRSLimit(const tensor3 & r,
		const tensor3 & gradient) const noexcept
{
	tensor3 UTCDFS;

	std::transform(r().cbegin(), r().cend(), UTCDFS.wr().begin(),
			[this](const auto r_j) 
			{	return this->UTCDFSLimiterCalculation(r_j);});

	return tensor3 { std::get<0>(UTCDFS()) * std::get<0>(gradient()),
			std::get<1>(UTCDFS()) * std::get<1>(gradient()), std::get<2>(
					UTCDFS()) * std::get<2>(gradient()), std::get<3>(UTCDFS())
					* std::get<3>(gradient()), std::get<4>(UTCDFS())
					* std::get<4>(gradient()), std::get<5>(UTCDFS())
					* std::get<5>(gradient()), std::get<6>(UTCDFS())
					* std::get<6>(gradient()), std::get<7>(UTCDFS())
					* std::get<7>(gradient()), std::get<8>(UTCDFS())
					* std::get<8>(gradient()), std::get<9>(UTCDFS())
					* std::get<9>(gradient()), std::get<10>(UTCDFS())
					* std::get<10>(gradient()), std::get<11>(UTCDFS())
					* std::get<11>(gradient()), std::get<12>(UTCDFS())
					* std::get<12>(gradient()), std::get<13>(UTCDFS())
					* std::get<13>(gradient()), std::get<14>(UTCDFS())
					* std::get<14>(gradient()), std::get<15>(UTCDFS())
					* std::get<15>(gradient()), std::get<16>(UTCDFS())
					* std::get<16>(gradient()), std::get<17>(UTCDFS())
					* std::get<17>(gradient()), std::get<18>(UTCDFS())
					* std::get<18>(gradient()), std::get<19>(UTCDFS())
					* std::get<19>(gradient()), std::get<20>(UTCDFS())
					* std::get<20>(gradient()), std::get<21>(UTCDFS())
					* std::get<21>(gradient()), std::get<22>(UTCDFS())
					* std::get<22>(gradient()), std::get<23>(UTCDFS())
					* std::get<23>(gradient()), std::get<24>(UTCDFS())
					* std::get<24>(gradient()), std::get<25>(UTCDFS())
					* std::get<25>(gradient()), std::get<26>(UTCDFS())
					* std::get<26>(gradient()) };
}
