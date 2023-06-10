/*
 * schemiTensor.hpp
 *
 *  Created on: 2019/11/12
 *      Author: Maxim Boldyrev
 *
 *      Class for 3D 2nd rank tensor.
 */

#ifndef SCHEMITENSOR_HPP_
#define SCHEMITENSOR_HPP_

#include <array>

#include "scalar.hpp"

namespace schemi
{
class tensor
{
public:
	tensor() noexcept = default;

	explicit tensor(const scalar inValue) noexcept;

	tensor(const scalar VXXin, const scalar VXYin, const scalar VXZin,
			const scalar VYXin, const scalar VYYin, const scalar VYZin,
			const scalar VZXin, const scalar VZYin, const scalar VZZin) noexcept;

	const std::array<scalar, 9>& v() const noexcept;

	std::array<scalar, 9>& v_r() noexcept;

	scalar trace() const noexcept;

	scalar mag() const noexcept;

	tensor& transpose() noexcept;

	tensor operator+(const tensor & inTensor) const noexcept;

	tensor& operator+=(const tensor & inTensor) noexcept;

	tensor operator-(const tensor & inTensor) const noexcept;

	tensor& operator-=(const tensor & inTensor) noexcept;

	tensor operator*(const scalar inScalar) const noexcept;

	tensor& operator*=(const scalar inScalar) noexcept;

	tensor operator/(const scalar inScalar) const noexcept;

	tensor& operator/=(const scalar inScalar) noexcept;

	tensor& operator=(const scalar inScalar) noexcept;

	constexpr static std::size_t vsize = 9;

private:
	std::array<scalar, 9> value =
			{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
};

tensor operator*(const scalar inScalar, const tensor & inTensor) noexcept;
}  // namespace schemi

#endif /* TENSOR_HPP_ */
