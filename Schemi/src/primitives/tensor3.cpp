/*
 * tensor3.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensor3.hpp"

#include <algorithm>
#include <cmath>

schemi::tensor3::tensor3(const scalar inValue) noexcept :
		value { inValue, inValue, inValue, inValue, inValue, inValue, inValue,
				inValue, inValue,

				inValue, inValue, inValue, inValue, inValue, inValue, inValue,
				inValue, inValue,

				inValue, inValue, inValue, inValue, inValue, inValue, inValue,
				inValue, inValue, }
{
}

schemi::tensor3::tensor3(const scalar VXXXin, const scalar VXXYin,
		const scalar VXXZin, const scalar VXYXin, const scalar VXYYin,
		const scalar VXYZin, const scalar VXZXin, const scalar VXZYin,
		const scalar VXZZin,

		const scalar VYXXin, const scalar VYXYin, const scalar VYXZin,
		const scalar VYYXin, const scalar VYYYin, const scalar VYYZin,
		const scalar VYZXin, const scalar VYZYin, const scalar VYZZin,

		const scalar VZXXin, const scalar VZXYin, const scalar VZXZin,
		const scalar VZYXin, const scalar VZYYin, const scalar VZYZin,
		const scalar VZZXin, const scalar VZZYin, const scalar VZZZin) noexcept :
		value { VXXXin, VXXYin, VXXZin, VXYXin, VXYYin, VXYZin, VXZXin, VXZYin,
				VXZZin,

				VYXXin, VYXYin, VYXZin, VYYXin, VYYYin, VYYZin, VYZXin, VYZYin,
				VYZZin,

				VZXXin, VZXYin, VZXZin, VZYXin, VZYYin, VZYZin, VZZXin, VZZYin,
				VZZZin }
{
}

const std::array<schemi::scalar, 27>& schemi::tensor3::operator()() const noexcept
{
	return value;
}

std::array<schemi::scalar, 27>& schemi::tensor3::r() noexcept
{
	return value;
}

schemi::scalar schemi::tensor3::trace() const noexcept
{
	return std::get<0>(value) + std::get<13>(value) + std::get<26>(value);
}

schemi::scalar schemi::tensor3::mag() const noexcept
{
	scalar result { 0 };

	for (std::size_t i = 0; i < value.size(); ++i)
		result += value[i] * value[i];

	return std::sqrt(result);
}

schemi::tensor3& schemi::tensor3::transpose() noexcept
{
	const auto buf(value);

	for (std::size_t i = 0; i < value.size(); ++i)
	{
		const auto d3 { i / 9 };
		const auto row { (i - 9 * d3) / 3 };
		const auto column { i - 9 * d3 - 3 * row };

		if (!((d3 == row) && (row == column)))
		{
			const std::size_t transposed_index = column * 9 + d3 * 3 + row;

			value[transposed_index] = buf[i];
		}
	}

	return *this;
}

schemi::tensor3 schemi::tensor3::operator+(
		const tensor3 & inTensor) const noexcept
{
	tensor3 sum(inTensor);
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.r()[i] += value[i];

	return sum;
}

schemi::tensor3& schemi::tensor3::operator+=(const tensor3 & inTensor) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] += inTensor()[i];

	return *this;
}

schemi::tensor3 schemi::tensor3::operator-(
		const tensor3 & inTensor) const noexcept
{
	tensor3 sum(inTensor * (-1.));
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.r()[i] += value[i];

	return sum;
}

schemi::tensor3& schemi::tensor3::operator-=(const tensor3 & inTensor) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] -= inTensor()[i];

	return *this;
}

schemi::tensor3 schemi::tensor3::operator*(const scalar inScalar) const noexcept
{
	return tensor3(std::get<0>(value) * inScalar, std::get<1>(value) * inScalar,
			std::get<2>(value) * inScalar, std::get<3>(value) * inScalar,
			std::get<4>(value) * inScalar, std::get<5>(value) * inScalar,
			std::get<6>(value) * inScalar, std::get<7>(value) * inScalar,
			std::get<8>(value) * inScalar, std::get<9>(value) * inScalar,
			std::get<10>(value) * inScalar, std::get<11>(value) * inScalar,
			std::get<12>(value) * inScalar, std::get<13>(value) * inScalar,
			std::get<14>(value) * inScalar, std::get<15>(value) * inScalar,
			std::get<16>(value) * inScalar, std::get<17>(value) * inScalar,
			std::get<18>(value) * inScalar, std::get<19>(value) * inScalar,
			std::get<20>(value) * inScalar, std::get<21>(value) * inScalar,
			std::get<22>(value) * inScalar, std::get<23>(value) * inScalar,
			std::get<24>(value) * inScalar, std::get<25>(value) * inScalar,
			std::get<26>(value) * inScalar);
}

schemi::tensor3& schemi::tensor3::operator*=(const scalar inScalar) noexcept
{
	std::transform(value.begin(), value.end(), value.begin(),
			[inScalar](const auto i) 
			{	return i *inScalar;});

	return *this;
}

schemi::tensor3 schemi::tensor3::operator/(const scalar inScalar) const noexcept
{
	return tensor3(std::get<0>(value) / inScalar, std::get<1>(value) / inScalar,
			std::get<2>(value) / inScalar, std::get<3>(value) / inScalar,
			std::get<4>(value) / inScalar, std::get<5>(value) / inScalar,
			std::get<6>(value) / inScalar, std::get<7>(value) / inScalar,
			std::get<8>(value) / inScalar, std::get<9>(value) / inScalar,
			std::get<10>(value) / inScalar, std::get<11>(value) / inScalar,
			std::get<12>(value) / inScalar, std::get<13>(value) / inScalar,
			std::get<14>(value) / inScalar, std::get<15>(value) / inScalar,
			std::get<16>(value) / inScalar, std::get<17>(value) / inScalar,
			std::get<18>(value) / inScalar, std::get<19>(value) / inScalar,
			std::get<20>(value) / inScalar, std::get<21>(value) / inScalar,
			std::get<22>(value) / inScalar, std::get<23>(value) / inScalar,
			std::get<24>(value) / inScalar, std::get<25>(value) / inScalar,
			std::get<26>(value) / inScalar);
}

schemi::tensor3& schemi::tensor3::operator/=(const scalar inScalar) noexcept
{
	std::transform(value.begin(), value.end(), value.begin(),
			[inScalar](const auto i) 
			{	return i /inScalar;});

	return *this;
}

schemi::tensor3& schemi::tensor3::operator=(const scalar inScalar) noexcept
{
	value.fill(inScalar);

	return *this;
}

schemi::tensor3 schemi::operator*(const scalar inScalar,
		const tensor3 & inTensor) noexcept
{
	return inTensor * inScalar;
}
