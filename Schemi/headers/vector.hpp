/*
 * schemiVector.hpp
 *
 *  Created on: 2019/11/12
 *      Author: Maxim Boldyrev
 *
 *      Class for 3D vector.
 */

#ifndef SCHEMIVECTOR_HPP_
#define SCHEMIVECTOR_HPP_

#include <array>

#include "scalar.hpp"

namespace schemi
{
class vector
{
public:
	vector() noexcept = default;

	explicit vector(const scalar inValue) noexcept;

	vector(const scalar VXin, const scalar VYin, const scalar VZin) noexcept;

	const std::array<scalar, 3>& operator()() const noexcept;

	std::array<scalar, 3>& r() noexcept;

	scalar mag() const noexcept;

	vector operator+(const vector & inVector) const noexcept;

	vector& operator+=(const vector & inVector) noexcept;

	vector operator-(const vector & inVector) const noexcept;

	vector& operator-=(const vector & inVector) noexcept;

	vector operator*(const scalar inScalar) const noexcept;

	vector& operator*=(const scalar inScalar) noexcept;

	vector operator/(const scalar inScalar) const noexcept;

	vector& operator/=(const scalar inScalar) noexcept;

	vector& operator=(const scalar inScalar) noexcept;

	constexpr static std::size_t vsize = 3;

private:
	std::array<scalar, 3> value = { 0.0, 0.0, 0.0 };
};

vector operator*(const scalar inScalar, const vector & inVector) noexcept;
}  // namespace schemi

#endif /* SCHEMIVECTOR_HPP_ */
