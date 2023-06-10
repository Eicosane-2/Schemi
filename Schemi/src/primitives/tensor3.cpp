/*
 * tensor3.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensor3.hpp"

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

const std::array<schemi::scalar, 27>& schemi::tensor3::v() const noexcept
{
	return value;
}

std::array<schemi::scalar, 27>& schemi::tensor3::v_r() noexcept
{
	return value;
}

schemi::scalar schemi::tensor3::trace() const noexcept
{
	return value[0] + value[13] + value[26];
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
		sum.v_r()[i] += value[i];

	return sum;
}

schemi::tensor3& schemi::tensor3::operator+=(const tensor3 & inTensor) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] += inTensor.v()[i];

	return *this;
}

schemi::tensor3 schemi::tensor3::operator-(
		const tensor3 & inTensor) const noexcept
{
	tensor3 sum(inTensor * (-1.));
	for (std::size_t i = 0; i < value.size(); ++i)
		sum.v_r()[i] += value[i];

	return sum;
}

schemi::tensor3& schemi::tensor3::operator-=(const tensor3 & inTensor) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] -= inTensor.v()[i];

	return *this;
}

schemi::tensor3 schemi::tensor3::operator*(const scalar inScalar) const noexcept
{
	return tensor3(value[0] * inScalar, value[1] * inScalar,
			value[2] * inScalar, value[3] * inScalar, value[4] * inScalar,
			value[5] * inScalar, value[6] * inScalar, value[7] * inScalar,
			value[8] * inScalar, value[9] * inScalar, value[10] * inScalar,
			value[11] * inScalar, value[12] * inScalar, value[13] * inScalar,
			value[14] * inScalar, value[15] * inScalar, value[16] * inScalar,
			value[17] * inScalar, value[18] * inScalar, value[19] * inScalar,
			value[20] * inScalar, value[21] * inScalar, value[22] * inScalar,
			value[23] * inScalar, value[24] * inScalar, value[25] * inScalar,
			value[26] * inScalar);
}

schemi::tensor3& schemi::tensor3::operator*=(const scalar inScalar) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] *= inScalar;

	return *this;
}

schemi::tensor3 schemi::tensor3::operator/(const scalar inScalar) const noexcept
{
	return tensor3(value[0] / inScalar, value[1] / inScalar,
			value[2] / inScalar, value[3] / inScalar, value[4] / inScalar,
			value[5] / inScalar, value[6] / inScalar, value[7] / inScalar,
			value[8] / inScalar, value[9] / inScalar, value[10] / inScalar,
			value[11] / inScalar, value[12] / inScalar, value[13] / inScalar,
			value[14] / inScalar, value[15] / inScalar, value[16] / inScalar,
			value[17] / inScalar, value[18] / inScalar, value[19] / inScalar,
			value[20] / inScalar, value[21] / inScalar, value[22] / inScalar,
			value[23] / inScalar, value[24] / inScalar, value[25] / inScalar,
			value[26] / inScalar);
}

schemi::tensor3& schemi::tensor3::operator/=(const scalar inScalar) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] /= inScalar;

	return *this;
}

schemi::tensor3& schemi::tensor3::operator=(const scalar inScalar) noexcept
{
	for (std::size_t i = 0; i < value.size(); ++i)
		value[i] = inScalar;

	return *this;
}

schemi::tensor3 schemi::operator*(const scalar inScalar,
		const tensor3 & inTensor) noexcept
{
	return inTensor * inScalar;
}
