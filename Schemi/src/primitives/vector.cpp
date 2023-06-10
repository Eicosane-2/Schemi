/*
 * vector.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vector.hpp"

#include <array>
#include <cmath>

schemi::vector::vector(const schemi::scalar inValue) noexcept :
		value { inValue, inValue, inValue }
{
}

schemi::vector::vector(const scalar VXin, const scalar VYin,
		const scalar VZin) noexcept :
		value { VXin, VYin, VZin }
{
}

const std::array<schemi::scalar, 3>& schemi::vector::v() const noexcept
{
	return value;
}

std::array<schemi::scalar, 3>& schemi::vector::v_r() noexcept
{
	return value;
}

schemi::scalar schemi::vector::mag() const noexcept
{
	//return std::sqrt(
	//		value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
	return std::hypot(value[0], value[1], value[2]); //XXX May be slow.
}

schemi::vector schemi::vector::operator+(const vector & inVector) const noexcept
{
	vector sum(inVector);
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.v_r()[i] += value[i];

	return sum;
}

schemi::vector& schemi::vector::operator+=(const vector & inVector) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] += inVector.v()[i];

	return *this;
}

schemi::vector schemi::vector::operator-(const vector & inVector) const noexcept
{
	vector sum(inVector * (-1.));
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.v_r()[i] += value[i];

	return sum;
}

schemi::vector& schemi::vector::operator-=(const vector & inVector) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] -= inVector.v()[i];

	return *this;
}

schemi::vector schemi::vector::operator*(const scalar inScalar) const noexcept
{
	return vector(value[0] * inScalar, value[1] * inScalar, value[2] * inScalar);
}

schemi::vector& schemi::vector::operator*=(const scalar inScalar) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] *= inScalar;

	return *this;
}

schemi::vector schemi::vector::operator/(const scalar inScalar) const noexcept
{
	return vector(value[0] / inScalar, value[1] / inScalar, value[2] / inScalar);
}

schemi::vector& schemi::vector::operator/=(const scalar inScalar) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] /= inScalar;

	return *this;
}

schemi::vector& schemi::vector::operator=(const scalar inScalar) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] = inScalar;

	return *this;
}

schemi::vector schemi::operator*(const scalar inScalar,
		const vector & inVector) noexcept
{
	return inVector * inScalar;
}
