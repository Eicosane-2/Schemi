/*
 * tensor.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensor.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "vector.hpp"
#include "tensor3.hpp"

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

std::array<schemi::scalar, 9>& schemi::tensor::wr() noexcept
{
	return value;
}

schemi::scalar schemi::tensor::trace() const noexcept
{
	return std::get<0>(value) + std::get<4>(value) + std::get<8>(value);
}

schemi::scalar schemi::tensor::mag() const noexcept
{
	const scalar result = std::accumulate(value.cbegin(), value.cend(), 0.0,
			[](const scalar null, const scalar i) 
			{
				return null + i*i;
			});

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
		sum.wr()[i] += value[i];

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
		sum.wr()[i] += value[i];

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
	std::transform(value.cbegin(), value.cend(), value.begin(),
			[inScalar](const auto i) 
			{	return i *inScalar;});

	return *this;
}

schemi::tensor3 schemi::tensor::operator*(
		const vector & inVector) const noexcept
{
	const auto & tens = *this;

	return tensor3 {

	std::get<0>(tens()) * std::get<0>(inVector()), std::get<0>(tens())
			* std::get<1>(inVector()), std::get<0>(tens())
			* std::get<2>(inVector()), std::get<1>(tens())
			* std::get<0>(inVector()), std::get<1>(tens())
			* std::get<1>(inVector()), std::get<1>(tens())
			* std::get<2>(inVector()), std::get<2>(tens())
			* std::get<0>(inVector()), std::get<2>(tens())
			* std::get<1>(inVector()), std::get<2>(tens())
			* std::get<2>(inVector()),

	std::get<3>(tens()) * std::get<0>(inVector()), std::get<3>(tens())
			* std::get<1>(inVector()), std::get<3>(tens())
			* std::get<2>(inVector()), std::get<4>(tens())
			* std::get<0>(inVector()), std::get<4>(tens())
			* std::get<1>(inVector()), std::get<4>(tens())
			* std::get<2>(inVector()), std::get<5>(tens())
			* std::get<0>(inVector()), std::get<5>(tens())
			* std::get<1>(inVector()), std::get<5>(tens())
			* std::get<2>(inVector()),

	std::get<6>(tens()) * std::get<0>(inVector()), std::get<6>(tens())
			* std::get<1>(inVector()), std::get<6>(tens())
			* std::get<2>(inVector()), std::get<7>(tens())
			* std::get<0>(inVector()), std::get<7>(tens())
			* std::get<1>(inVector()), std::get<7>(tens())
			* std::get<2>(inVector()), std::get<8>(tens())
			* std::get<0>(inVector()), std::get<8>(tens())
			* std::get<1>(inVector()), std::get<8>(tens())
			* std::get<2>(inVector())

	};
}

schemi::vector schemi::tensor::operator&(const vector & inVector) const noexcept
{
	const auto & tens = *this;

	return vector { std::get<0>(tens()) * std::get<0>(inVector())
			+ std::get<1>(tens()) * std::get<1>(inVector())
			+ std::get<2>(tens()) * std::get<2>(inVector()), std::get<3>(tens())
			* std::get<0>(inVector())
			+ std::get<4>(tens()) * std::get<1>(inVector())
			+ std::get<5>(tens()) * std::get<2>(inVector()), std::get<6>(tens())
			* std::get<0>(inVector())
			+ std::get<7>(tens()) * std::get<1>(inVector())
			+ std::get<8>(tens()) * std::get<2>(inVector()) };
}

schemi::tensor schemi::tensor::operator&(const tensor & inTensor) const noexcept
{
	const auto & tens = *this;

	return tensor { std::get<0>(tens()) * std::get<0>(inTensor())
			+ std::get<1>(tens()) * std::get<3>(inTensor())
			+ std::get<2>(tens()) * std::get<6>(inTensor()), std::get<0>(tens())
			* std::get<1>(inTensor())
			+ std::get<1>(tens()) * std::get<4>(inTensor())
			+ std::get<2>(tens()) * std::get<7>(inTensor()), std::get<0>(tens())
			* std::get<2>(inTensor())
			+ std::get<1>(tens()) * std::get<5>(inTensor())
			+ std::get<2>(tens()) * std::get<8>(inTensor()), std::get<3>(tens())
			* std::get<0>(inTensor())
			+ std::get<4>(tens()) * std::get<3>(inTensor())
			+ std::get<5>(tens()) * std::get<6>(inTensor()), std::get<3>(tens())
			* std::get<1>(inTensor())
			+ std::get<4>(tens()) * std::get<4>(inTensor())
			+ std::get<5>(tens()) * std::get<7>(inTensor()), std::get<3>(tens())
			* std::get<2>(inTensor())
			+ std::get<4>(tens()) * std::get<5>(inTensor())
			+ std::get<5>(tens()) * std::get<8>(inTensor()), std::get<6>(tens())
			* std::get<0>(inTensor())
			+ std::get<7>(tens()) * std::get<3>(inTensor())
			+ std::get<8>(tens()) * std::get<6>(inTensor()), std::get<6>(tens())
			* std::get<1>(inTensor())
			+ std::get<7>(tens()) * std::get<4>(inTensor())
			+ std::get<8>(tens()) * std::get<7>(inTensor()), std::get<6>(tens())
			* std::get<2>(inTensor())
			+ std::get<7>(tens()) * std::get<5>(inTensor())
			+ std::get<8>(tens()) * std::get<8>(inTensor()) };
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
	std::transform(value.cbegin(), value.cend(), value.begin(),
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

schemi::vector schemi::operator&(const vector & inVector,
		const tensor & inTensor) noexcept
{
	tensor bufTensor(inTensor);
	bufTensor.transpose();

	return bufTensor & inVector;
}

schemi::tensor3 schemi::operator*(const vector & inVector,
		const tensor & inTensor) noexcept
{
	return tensor3 {

	std::get<0>(inVector()) * std::get<0>(inTensor()), std::get<0>(inVector())
			* std::get<1>(inTensor()), std::get<0>(inVector())
			* std::get<2>(inTensor()), std::get<0>(inVector())
			* std::get<3>(inTensor()), std::get<0>(inVector())
			* std::get<4>(inTensor()), std::get<0>(inVector())
			* std::get<5>(inTensor()), std::get<0>(inVector())
			* std::get<6>(inTensor()), std::get<0>(inVector())
			* std::get<7>(inTensor()), std::get<0>(inVector())
			* std::get<8>(inTensor()),

	std::get<1>(inVector()) * std::get<0>(inTensor()), std::get<1>(inVector())
			* std::get<1>(inTensor()), std::get<1>(inVector())
			* std::get<2>(inTensor()), std::get<1>(inVector())
			* std::get<3>(inTensor()), std::get<1>(inVector())
			* std::get<4>(inTensor()), std::get<1>(inVector())
			* std::get<5>(inTensor()), std::get<1>(inVector())
			* std::get<6>(inTensor()), std::get<1>(inVector())
			* std::get<7>(inTensor()), std::get<1>(inVector())
			* std::get<8>(inTensor()),

	std::get<2>(inVector()) * std::get<0>(inTensor()), std::get<2>(inVector())
			* std::get<1>(inTensor()), std::get<2>(inVector())
			* std::get<2>(inTensor()), std::get<2>(inVector())
			* std::get<3>(inTensor()), std::get<2>(inVector())
			* std::get<4>(inTensor()), std::get<2>(inVector())
			* std::get<5>(inTensor()), std::get<2>(inVector())
			* std::get<6>(inTensor()), std::get<2>(inVector())
			* std::get<7>(inTensor()), std::get<2>(inVector())
			* std::get<8>(inTensor())

	};
}
