/*
 * tensor.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensor.hpp"

#include <algorithm>
#include <cmath>

schemi::tensor::tensor(const scalar inValue) noexcept :
		value { inValue, inValue, inValue, inValue, inValue, inValue, inValue,
				inValue, inValue }
{
}

schemi::tensor::tensor(const scalar VXXin, const scalar VXYin,
		const scalar VXZin, const scalar VYXin, const scalar VYYin,
		const scalar VYZin, const scalar VZXin, const scalar VZYin,
		const scalar VZZin) noexcept :
		value { VXXin, VXYin, VXZin, VYXin, VYYin, VYZin, VZXin, VZYin, VZZin }
{
}

const std::array<schemi::scalar, 9>& schemi::tensor::operator()() const noexcept
{
	return value;
}

std::array<schemi::scalar, 9>& schemi::tensor::r() noexcept
{
	return value;
}

schemi::scalar schemi::tensor::trace() const noexcept
{
	return std::get<0>(value) + std::get<4>(value) + std::get<8>(value);
}

schemi::scalar schemi::tensor::mag() const noexcept
{
	scalar result { 0 };

	for (std::size_t i = 0; i < value.size(); ++i)
		result += value[i] * value[i];

	return std::sqrt(result);
}

schemi::tensor& schemi::tensor::transpose() noexcept
{
	const auto buf(value);

	std::get<1>(value) = std::get<3>(buf);
	std::get<2>(value) = std::get<6>(buf);
	std::get<3>(value) = std::get<1>(buf);
	std::get<5>(value) = std::get<7>(buf);
	std::get<6>(value) = std::get<2>(buf);
	std::get<7>(value) = std::get<5>(buf);

	return *this;
}

schemi::tensor schemi::tensor::operator+(const tensor & inTensor) const noexcept
{
	tensor sum(inTensor);
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.r()[i] += value[i];

	return sum;
}

schemi::tensor& schemi::tensor::operator+=(const tensor & inTensor) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] += inTensor()[i];

	return *this;
}

schemi::tensor schemi::tensor::operator-(const tensor & inTensor) const noexcept
{
	tensor sum(inTensor * (-1.));
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.r()[i] += value[i];

	return sum;
}

schemi::tensor& schemi::tensor::operator-=(const tensor & inTensor) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] -= inTensor()[i];

	return *this;
}

schemi::tensor schemi::tensor::operator*(const scalar inScalar) const noexcept
{
	return tensor(std::get<0>(value) * inScalar, std::get<1>(value) * inScalar,
			std::get<2>(value) * inScalar, std::get<3>(value) * inScalar,
			std::get<4>(value) * inScalar, std::get<5>(value) * inScalar,
			std::get<6>(value) * inScalar, std::get<7>(value) * inScalar,
			std::get<8>(value) * inScalar);
}

schemi::tensor& schemi::tensor::operator*=(const scalar inScalar) noexcept
{
	std::transform(value.begin(), value.end(), value.begin(),
			[inScalar](const auto i) 
			{	return i *inScalar;});

	return *this;
}

schemi::tensor schemi::tensor::operator/(const scalar inScalar) const noexcept
{
	return tensor(std::get<0>(value) / inScalar, std::get<1>(value) / inScalar,
			std::get<2>(value) / inScalar, std::get<3>(value) / inScalar,
			std::get<4>(value) / inScalar, std::get<5>(value) / inScalar,
			std::get<6>(value) / inScalar, std::get<7>(value) / inScalar,
			std::get<8>(value) / inScalar);
}

schemi::tensor& schemi::tensor::operator/=(const scalar inScalar) noexcept
{
	std::transform(value.begin(), value.end(), value.begin(),
			[inScalar](const auto i) 
			{	return i /inScalar;});

	return *this;
}

schemi::tensor& schemi::tensor::operator=(const scalar inScalar) noexcept
{
	value.fill(inScalar);

	return *this;
}

schemi::tensor schemi::operator*(const scalar inScalar,
		const tensor & inTensor) noexcept
{
	return inTensor * inScalar;
}
