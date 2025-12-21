/*
 * fieldOperations.hpp
 *
 *  Created on: 2025/11/22
 *      Author: Maxim Boldyrev
 */

#ifndef FIELDOPERATIONS_HPP_
#define FIELDOPERATIONS_HPP_

#include <algorithm>

#include "field.hpp"
#include "mathObjects.hpp"
#include "doubleDotProduct.hpp"

namespace schemi
{
/*+*/
template<mathObjects T, typename typeOfField>
auto operator+(const field<T, typeOfField> & F1,
		const field<T, typeOfField> & F2) noexcept
{
	field<T, typeOfField> retField
	{	F1.meshRef(), T
		{	0}};

	retField.val() = F1.cval() + F2.cval();

	return retField;
}

template<mathObjects T, typename typeOfField>
auto operator+(const field<T, typeOfField> & F1, const T & F2) noexcept
{
	field<T, typeOfField> retField
	{	F1.meshRef(), T
		{	0}};

	retField.val() = F1.cval() + F2;

	return F1;
}

template<mathObjects T, typename typeOfField>
auto operator+(const T & F2, const field<T, typeOfField> & F1) noexcept
{
	field<T, typeOfField> retField
	{	F1.meshRef(), T
		{	0}};

	retField.val() = F2 + F1.cval();

	return F1;
}

template<mathObjects T, typename typeOfField>
auto& operator+=(field<T, typeOfField> & F1,
		const field<T, typeOfField> & F2) noexcept
{
	F1.val() += F2.cval();

	return F1;
}

template<mathObjects T, typename typeOfField>
auto& operator+=(field<T, typeOfField> & F1, const T & F2) noexcept
{
	F1.val() += F2;

	return F1;
}
/*+*/

/*-*/
template<mathObjects T, typename typeOfField>
auto operator-(const field<T, typeOfField> & F1,
		const field<T, typeOfField> & F2) noexcept
{
	field<T, typeOfField> retField
	{	F1.meshRef(), T
		{	0}};

	retField.val() = F1.cval() - F2.cval();

	return retField;
}

template<mathObjects T, typename typeOfField>
auto operator-(const field<T, typeOfField> & F1, const T & F2) noexcept
{
	field<T, typeOfField> retField
	{	F1.meshRef(), T
		{	0}};

	retField.val() = F1.cval() - F2;

	return retField;
}

template<mathObjects T, typename typeOfField>
auto operator-(const T & F2, const field<T, typeOfField> & F1) noexcept
{
	field<T, typeOfField> retField
	{	F1.meshRef(), T
		{	0}};

	retField.val() = F2 - F1.cval();

	return retField;
}

template<mathObjects T, typename typeOfField>
auto& operator-=(field<T, typeOfField> & F1,
		const field<T, typeOfField> & F2) noexcept
{
	F1.val() -= F2.cval();

	return F1;
}

template<mathObjects T, typename typeOfField>
auto& operator-=(field<T, typeOfField> & F1, const T & F2) noexcept
{
	F1.val() -= F2;

	return F1;
}
/*-*/

/***/
template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator*(const field<T1, typeOfField> & F1,
		const field<T2, typeOfField> & F2) noexcept
{
	using resultType = decltype(T1() * T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F2.cval()), std::begin(retField.val()),
			[](const auto & F1_i, const auto & F2_i) -> resultType
			{
				return F1_i * F2_i;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator*(const field<T1, typeOfField> & F1, const T2 & F2) noexcept
{
	using resultType = decltype(T1() * T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F1_i * F2;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator*(const T2 & F2, const field<T1, typeOfField> & F1) noexcept
{
	using resultType = decltype(T2() * T1());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F2 * F1_i;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto& operator*=(field<T1, typeOfField> & F1,
		const field<T2, typeOfField> & F2) noexcept
{
	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F2.cval()), std::begin(F1.val()),
			[](const auto & F1_i, const auto & F2_i) -> T1
			{
				return F1_i * F2_i;
			});

	return F1;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto& operator*=(field<T1, typeOfField> & F1, const T2 & F2) noexcept
{
	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F1.val()), [&F2](const auto & F1_i) -> T1
			{
				return F1_i * F2;
			});

	return F1;
}
/***/

/*/*/
template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator/(const field<T1, typeOfField> & F1,
		const field<T2, typeOfField> & F2) noexcept
{
	using resultType = decltype(T1() / T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F2.cval()), std::begin(retField.val()),
			[](const auto & F1_i, const auto & F2_i) -> resultType
			{
				return F1_i / F2_i;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator/(const field<T1, typeOfField> & F1, const T2 & F2) noexcept
{
	using resultType = decltype(T1() / T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F1_i / F2;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator/(const T2 & F2, const field<T1, typeOfField> & F1) noexcept
{
	using resultType = decltype(T2() / T1());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F2 / F1_i;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto& operator/=(field<T1, typeOfField> & F1,
		field<T2, typeOfField> & F2) noexcept
{
	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F2.cval()), std::begin(F1.val()),
			[](const auto & F1_i, const auto & F2_i) -> T1
			{
				return F1_i / F2_i;
			});

	return F1;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto& operator/=(field<T1, typeOfField> & F1, const T2 & F2) noexcept
{
	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F1.val()), [&F2](const auto & F1_i) -> T1
			{
				return F1_i / F2;
			});

	return F1;
}
/*/*/

/*&*/
template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator&(const field<T1, typeOfField> & F1,
		const field<T2, typeOfField> & F2) noexcept
{
	using resultType = decltype(T1() & T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F2.cval()), std::begin(retField.val()),
			[](const auto & F1_i, const auto & F2_i) -> resultType
			{
				return F1_i & F2_i;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator&(const field<T1, typeOfField> & F1, const T2 & F2) noexcept
{
	using resultType = decltype(T1() & T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F1_i & F2;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator&(const T2 & F2, const field<T1, typeOfField> & F1) noexcept
{
	using resultType = decltype(T2() & T1());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F2 & F1_i;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto& operator&=(field<T1, typeOfField> & F1,
		const field<T2, typeOfField> & F2) noexcept
{
	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F2.cval()), std::begin(F1.val()),
			[](const auto & F1_i, const auto & F2_i) -> T1
			{
				return F1_i & F2_i;
			});

	return F1;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto& operator&=(field<T1, typeOfField> & F1, const T2 & F2) noexcept
{
	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F1.val()), [&F2](const auto & F1_i) -> T1
			{
				return F1_i & F2;
			});

	return F1;
}
/*&*/

/*&&*/
template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator&&(const field<T1, typeOfField> & F1,
		const field<T2, typeOfField> & F2) noexcept
{
	using resultType = decltype(T1() && T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F2.cval()), std::begin(retField.val()),
			[](const auto & F1_i, const auto & F2_i) -> resultType
			{
				return F1_i && F2_i;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator&&(const field<T1, typeOfField> & F1, const T2 & F2) noexcept
{
	using resultType = decltype(T1() && T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F1_i && F2;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator&&(const T2 & F2, const field<T1, typeOfField> & F1) noexcept
{
	using resultType = decltype(T2() && T1());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F2 && F1_i;
			});

	return retField;
}
/*&&*/

/*^*/
template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator^(const field<T1, typeOfField> & F1,
		const field<T2, typeOfField> & F2) noexcept
{
	using resultType = decltype(T1() ^ T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F2.cval()), std::begin(retField.val()),
			[](const auto & F1_i, const auto & F2_i) -> resultType
			{
				return F1_i ^ F2_i;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator^(const field<T1, typeOfField> & F1, const T2 & F2) noexcept
{
	using resultType = decltype(T1() ^ T2());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F1_i ^ F2;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto operator^(const T2 & F2, const field<T1, typeOfField> & F1) noexcept
{
	using resultType = decltype(T2() ^ T1());

	field<resultType, typeOfField> retField
	{	F1.meshRef(), resultType
		{	0}};

	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(retField.val()), [&F2](const auto & F1_i) -> resultType
			{
				return F2 ^ F1_i;
			});

	return retField;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto& operator^=(const field<T1, typeOfField> & F1,
		const field<T2, typeOfField> & F2) noexcept
{
	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F2.cval()), std::begin(F1.val()),
			[](const auto & F1_i, const auto & F2_i) -> T1
			{
				return F1_i ^ F2_i;
			});

	return F1;
}

template<mathObjects T1, mathObjects T2, typename typeOfField>
auto& operator^=(const field<T1, typeOfField> & F1, const T2 & F2) noexcept
{
	std::transform(std::begin(F1.cval()), std::end(F1.cval()),
			std::begin(F1.val()), [&F2](const auto & F1_i) -> T1
			{
				return F1_i ^ F2;
			});

	return F1;
}
/*^*/
}

#endif /* FIELDOPERATIONS_HPP_ */
