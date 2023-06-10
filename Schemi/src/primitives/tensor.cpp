/*
 * tensor.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensor.hpp"

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

const std::array<schemi::scalar, 9>& schemi::tensor::v() const noexcept
{
	return value;
}

std::array<schemi::scalar, 9>& schemi::tensor::v_r() noexcept
{
	return value;
}

schemi::scalar schemi::tensor::trace() const noexcept
{
	return value[0] + value[4] + value[8];
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

	value[1] = buf[3];
	value[2] = buf[6];
	value[3] = buf[1];
	value[5] = buf[7];
	value[6] = buf[2];
	value[7] = buf[5];

	return *this;
}

schemi::tensor schemi::tensor::operator+(const tensor & inTensor) const noexcept
{
	tensor sum(inTensor);
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.v_r()[i] += value[i];

	return sum;
}

schemi::tensor& schemi::tensor::operator+=(const tensor & inTensor) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] += inTensor.v()[i];

	return *this;
}

schemi::tensor schemi::tensor::operator-(const tensor & inTensor) const noexcept
{
	tensor sum(inTensor * (-1.));
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.v_r()[i] += value[i];

	return sum;
}

schemi::tensor& schemi::tensor::operator-=(const tensor & inTensor) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] -= inTensor.v()[i];

	return *this;
}

schemi::tensor schemi::tensor::operator*(const scalar inScalar) const noexcept
{
	return tensor(value[0] * inScalar, value[1] * inScalar, value[2] * inScalar,
			value[3] * inScalar, value[4] * inScalar, value[5] * inScalar,
			value[6] * inScalar, value[7] * inScalar, value[8] * inScalar);
}

schemi::tensor& schemi::tensor::operator*=(const scalar inScalar) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] *= inScalar;

	return *this;
}

schemi::tensor schemi::tensor::operator/(const scalar inScalar) const noexcept
{
	return tensor(value[0] / inScalar, value[1] / inScalar, value[2] / inScalar,
			value[3] / inScalar, value[4] / inScalar, value[5] / inScalar,
			value[6] / inScalar, value[7] / inScalar, value[8] / inScalar);
}

schemi::tensor& schemi::tensor::operator/=(const scalar inScalar) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] /= inScalar;

	return *this;
}

schemi::tensor& schemi::tensor::operator=(const scalar inScalar) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] = inScalar;

	return *this;
}

schemi::tensor schemi::operator*(const scalar inScalar,
		const tensor & inTensor) noexcept
{
	return inTensor * inScalar;
}
