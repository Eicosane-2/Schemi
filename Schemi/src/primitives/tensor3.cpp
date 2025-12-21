/*
 * tensor3.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensor3.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "vector.hpp"
#include "tensor.hpp"

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

std::array<schemi::scalar, 27>& schemi::tensor3::wr() noexcept
{
	return value;
}

schemi::scalar schemi::tensor3::trace() const noexcept
{
	return std::get<0>(value) + std::get<13>(value) + std::get<26>(value);
}

schemi::scalar schemi::tensor3::mag() const noexcept
{
	const scalar result = std::accumulate(value.cbegin(), value.cend(), 0.0,
			[](const scalar null, const scalar i) 
			{
				return null + i*i;
			});

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
		sum.wr()[i] += value[i];

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
		sum.wr()[i] += value[i];

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
	std::transform(value.cbegin(), value.cend(), value.begin(),
			[inScalar](const auto i) 
			{	return i *inScalar;});

	return *this;
}

schemi::tensor schemi::tensor3::operator&(
		const vector & inVector) const noexcept
{
	const auto tens3 = *this;

	return tensor { std::get<0>(tens3()) * std::get<0>(inVector())
			+ std::get<1>(tens3()) * std::get<1>(inVector())
			+ std::get<2>(tens3()) * std::get<2>(inVector()), std::get<3>(
			tens3()) * std::get<0>(inVector())
			+ std::get<4>(tens3()) * std::get<1>(inVector())
			+ std::get<5>(tens3()) * std::get<2>(inVector()), std::get<6>(
			tens3()) * std::get<0>(inVector())
			+ std::get<7>(tens3()) * std::get<1>(inVector())
			+ std::get<8>(tens3()) * std::get<2>(inVector()),

	std::get<9 + 0>(tens3()) * std::get<0>(inVector())
			+ std::get<9 + 1>(tens3()) * std::get<1>(inVector())
			+ std::get<9 + 2>(tens3()) * std::get<2>(inVector()),
			std::get<9 + 3>(tens3()) * std::get<0>(inVector())
					+ std::get<9 + 4>(tens3()) * std::get<1>(inVector())
					+ std::get<9 + 5>(tens3()) * std::get<2>(inVector()),
			std::get<9 + 6>(tens3()) * std::get<0>(inVector())
					+ std::get<9 + 7>(tens3()) * std::get<1>(inVector())
					+ std::get<9 + 8>(tens3()) * std::get<2>(inVector()),

			std::get<18 + 0>(tens3()) * std::get<0>(inVector())
					+ std::get<18 + 1>(tens3()) * std::get<1>(inVector())
					+ std::get<18 + 2>(tens3()) * std::get<2>(inVector()),
			std::get<18 + 3>(tens3()) * std::get<0>(inVector())
					+ std::get<18 + 4>(tens3()) * std::get<1>(inVector())
					+ std::get<18 + 5>(tens3()) * std::get<2>(inVector()),
			std::get<18 + 6>(tens3()) * std::get<0>(inVector())
					+ std::get<18 + 7>(tens3()) * std::get<1>(inVector())
					+ std::get<18 + 8>(tens3()) * std::get<2>(inVector()) };
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
	std::transform(value.cbegin(), value.cend(), value.begin(),
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

schemi::tensor schemi::operator&(const vector & inVector,
		const tensor3 & inTensor) noexcept
{
	return tensor { std::get<0>(inTensor()) * std::get<0>(inVector())
			+ std::get<3>(inTensor()) * std::get<1>(inVector())
			+ std::get<6>(inTensor()) * std::get<2>(inVector()), std::get<1>(
			inTensor()) * std::get<0>(inVector())
			+ std::get<4>(inTensor()) * std::get<1>(inVector())
			+ std::get<7>(inTensor()) * std::get<2>(inVector()), std::get<2>(
			inTensor()) * std::get<0>(inVector())
			+ std::get<5>(inTensor()) * std::get<1>(inVector())
			+ std::get<8>(inTensor()) * std::get<2>(inVector()),

	std::get<9 + 0>(inTensor()) * std::get<0>(inVector())
			+ std::get<9 + 3>(inTensor()) * std::get<1>(inVector())
			+ std::get<9 + 6>(inTensor()) * std::get<2>(inVector()), std::get<
			9 + 1>(inTensor()) * std::get<0>(inVector())
			+ std::get<9 + 4>(inTensor()) * std::get<1>(inVector())
			+ std::get<9 + 7>(inTensor()) * std::get<2>(inVector()), std::get<
			9 + 2>(inTensor()) * std::get<0>(inVector())
			+ std::get<9 + 5>(inTensor()) * std::get<1>(inVector())
			+ std::get<9 + 8>(inTensor()) * std::get<2>(inVector()),

	std::get<18 + 0>(inTensor()) * std::get<0>(inVector())
			+ std::get<18 + 3>(inTensor()) * std::get<1>(inVector())
			+ std::get<18 + 6>(inTensor()) * std::get<2>(inVector()), std::get<
			18 + 1>(inTensor()) * std::get<0>(inVector())
			+ std::get<18 + 4>(inTensor()) * std::get<1>(inVector())
			+ std::get<18 + 7>(inTensor()) * std::get<2>(inVector()), std::get<
			18 + 2>(inTensor()) * std::get<0>(inVector())
			+ std::get<18 + 5>(inTensor()) * std::get<1>(inVector())
			+ std::get<18 + 8>(inTensor()) * std::get<2>(inVector()) };
}
