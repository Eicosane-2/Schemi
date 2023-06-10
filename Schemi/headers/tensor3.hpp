/*
 * tensor3.hpp
 *
 *  Created on: 2021/05/27
 *      Author: Maxim Boldyrev
 *
 *      Class for 3D 3d rank tensor.
 */

#ifndef TENSOR3_HPP_
#define TENSOR3_HPP_

#include <array>

#include "scalar.hpp"

namespace schemi
{
class tensor3
{
public:
	tensor3() noexcept = default;

	explicit tensor3(const scalar inValue) noexcept;

	tensor3(const scalar VXXXin, const scalar VXXYin, const scalar VXXZin,
			const scalar VXYXin, const scalar VXYYin, const scalar VXYZin,
			const scalar VXZXin, const scalar VXZYin, const scalar VXZZin,

			const scalar VYXXin, const scalar VYXYin, const scalar VYXZin,
			const scalar VYYXin, const scalar VYYYin, const scalar VYYZin,
			const scalar VYZXin, const scalar VYZYin, const scalar VYZZin,

			const scalar VZXXin, const scalar VZXYin, const scalar VZXZin,
			const scalar VZYXin, const scalar VZYYin, const scalar VZYZin,
			const scalar VZZXin, const scalar VZZYin,
			const scalar VZZZin) noexcept;

	const std::array<scalar, 27>& v() const noexcept;

	std::array<scalar, 27>& v_r() noexcept;

	scalar trace() const noexcept;

	scalar mag() const noexcept;

	tensor3& transpose() noexcept;

	tensor3 operator+(const tensor3 & inTensor) const noexcept;

	tensor3& operator+=(const tensor3 & inTensor) noexcept;

	tensor3 operator-(const tensor3 & inTensor) const noexcept;

	tensor3& operator-=(const tensor3 & inTensor) noexcept;

	tensor3 operator*(const scalar inScalar) const noexcept;

	tensor3& operator*=(const scalar inScalar) noexcept;

	tensor3 operator/(const scalar inScalar) const noexcept;

	tensor3& operator/=(const scalar inScalar) noexcept;

	tensor3& operator=(const scalar inScalar) noexcept;

	constexpr static std::size_t vsize = 27;

private:
	std::array<scalar, 27> value = {

	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
};

tensor3 operator*(const scalar inScalar, const tensor3 & inTensor) noexcept;
}  // namespace schemi

#endif /* TENSOR3_HPP_ */
