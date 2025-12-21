/*
 * vector.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "tensor.hpp"

schemi::vector::vector(const schemi::scalar inValue) noexcept :
		value { inValue, inValue, inValue }
{
}

schemi::vector::vector(const scalar VXin, const scalar VYin,
		const scalar VZin) noexcept :
		value { VXin, VYin, VZin }
{
}

const std::array<schemi::scalar, 3>& schemi::vector::operator()() const noexcept
{
	return value;
}

std::array<schemi::scalar, 3>& schemi::vector::wr() noexcept
{
	return value;
}

schemi::scalar schemi::vector::mag() const noexcept
{
	//return std::sqrt(
	//		value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
	return std::hypot(std::get<0>(value), std::get<1>(value),
			std::get<2>(value)); //XXX May be slow.
}

schemi::vector schemi::vector::operator+(const vector & inVector) const noexcept
{
	vector sum(inVector);
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.wr()[i] += value[i];

	return sum;
}

schemi::vector& schemi::vector::operator+=(const vector & inVector) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] += inVector()[i];

	return *this;
}

schemi::vector schemi::vector::operator-(const vector & inVector) const noexcept
{
	vector sum(inVector * (-1.));
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.wr()[i] += value[i];

	return sum;
}

schemi::vector& schemi::vector::operator-=(const vector & inVector) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] -= inVector()[i];

	return *this;
}

schemi::vector schemi::vector::operator*(const scalar inScalar) const noexcept
{
	return vector(std::get<0>(value) * inScalar, std::get<1>(value) * inScalar,
			std::get<2>(value) * inScalar);
}

schemi::vector& schemi::vector::operator*=(const scalar inScalar) noexcept
{
	std::transform(value.cbegin(), value.cend(), value.begin(),
			[inScalar](const auto i) 
			{	return i *inScalar;});

	return *this;
}

schemi::tensor schemi::vector::operator*(const vector & inVector) const noexcept
{
	const auto & vec = *this;

	return tensor { std::get<0>(vec()) * std::get<0>(inVector()), std::get<0>(
			vec()) * std::get<1>(inVector()), std::get<0>(vec())
			* std::get<2>(inVector()), std::get<1>(vec())
			* std::get<0>(inVector()), std::get<1>(vec())
			* std::get<1>(inVector()), std::get<1>(vec())
			* std::get<2>(inVector()), std::get<2>(vec())
			* std::get<0>(inVector()), std::get<2>(vec())
			* std::get<1>(inVector()), std::get<2>(vec())
			* std::get<2>(inVector()) };
}

schemi::scalar schemi::vector::operator&(const vector & inVector) const noexcept
{
	const auto & vec = *this;

	scalar result = 0;

	for (std::size_t i = 0; i < vector::vsize; ++i)
		result += vec()[i] * inVector()[i];

	return result;
}

schemi::vector schemi::vector::operator^(const vector & inVector) const noexcept
{
	const auto & a = *this;
	const auto & b = inVector;

	return vector { std::get<1>(a()) * std::get<2>(b())
			- std::get<2>(a()) * std::get<1>(b()), -(std::get<0>(a())
			* std::get<2>(b()) - std::get<2>(a()) * std::get<0>(b())), std::get<
			0>(a()) * std::get<1>(b()) - std::get<1>(a()) * std::get<0>(b()) };
}

schemi::vector& schemi::vector::operator^=(const vector & inVector) noexcept
{
	const auto & a = *this;
	const auto & b = inVector;

	const auto result = vector { std::get<1>(a()) * std::get<2>(b())
			- std::get<2>(a()) * std::get<1>(b()), -(std::get<0>(a())
			* std::get<2>(b()) - std::get<2>(a()) * std::get<0>(b())), std::get<
			0>(a()) * std::get<1>(b()) - std::get<1>(a()) * std::get<0>(b()) };

	value = result();

	return *this;
}

schemi::vector schemi::vector::operator/(const scalar inScalar) const noexcept
{
	return vector(std::get<0>(value) / inScalar, std::get<1>(value) / inScalar,
			std::get<2>(value) / inScalar);
}

schemi::vector& schemi::vector::operator/=(const scalar inScalar) noexcept
{
	std::transform(value.cbegin(), value.cend(), value.begin(),
			[inScalar](const auto i) 
			{	return i /inScalar;});

	return *this;
}

schemi::vector& schemi::vector::operator=(const scalar inScalar) noexcept
{
	value.fill(inScalar);

	return *this;
}

schemi::vector schemi::operator*(const scalar inScalar,
		const vector & inVector) noexcept
{
	return inVector * inScalar;
}
