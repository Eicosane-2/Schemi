/*
 * schemiFieldProducts.hpp
 *
 *  Created on: 2021/02/10
 *      Author: Maxim Boldyrev
 *
 *      Functions concealing fields' products: "*", "&", "&&", "/".
 */

#ifndef FIELDPRODUCTS_HPP_
#define FIELDPRODUCTS_HPP_

#include "doubleDotProduct.hpp"
#include "dyadicProduct.hpp"
#include "tensorVectorDotProduct.hpp"
#include "vectorVectorDotProduct.hpp"

namespace schemi
{
template<typename arg1, typename arg2, typename typeOfField>
auto astProduct(const field<arg1, typeOfField> & F1,
		const field<arg2, typeOfField> & F2) noexcept
{
	using resultType = decltype(arg1() * arg2());

	field<resultType, typeOfField> retField { F1.meshRef(), resultType { 0 } };

	for (std::size_t i = 0; i < F1.size(); ++i)
		retField.ref_r()[i] = F1.ref()[i] * F2.ref()[i];

	return retField;
}

template<typename arg1, typename arg2, typename typeOfField>
auto astProduct(const field<arg1, typeOfField> & F1, const arg2 & F2) noexcept
{
	using resultType = decltype(arg1() * F2);

	field<resultType, typeOfField> retField { F1.meshRef(), resultType { 0 } };

	for (std::size_t i = 0; i < F1.size(); ++i)
		retField.ref_r()[i] = F1.ref()[i] * F2;

	return retField;
}

template<typename arg1, typename arg2, typename typeOfField>
auto& astProductSelf(field<arg1, typeOfField> & F1,
		const field<arg2, typeOfField> & F2) noexcept
{
	for (std::size_t i = 0; i < F1.size(); ++i)
		F1.ref_r()[i] *= F2.ref()[i];

	return F1;
}

template<typename arg1, typename arg2, typename typeOfField>
auto& astProductSelf(field<arg1, typeOfField> & F1, const arg2 & F2) noexcept
{
	for (std::size_t i = 0; i < F1.size(); ++i)
		F1.ref_r()[i] *= F2;

	return F1;
}

template<typename arg1, typename arg2, typename typeOfField>
auto ampProduct(const field<arg1, typeOfField> & F1,
		const field<arg2, typeOfField> & F2) noexcept
{
	using resultType = decltype(arg1() & arg2());

	field<resultType, typeOfField> retField { F1.meshRef(), resultType { 0 } };

	for (std::size_t i = 0; i < F1.size(); ++i)
		retField.ref_r()[i] = F1.ref()[i] & F2.ref()[i];

	return retField;
}

template<typename arg1, typename arg2, typename typeOfField>
auto ampProduct(const field<arg1, typeOfField> & F1, const arg2 & F2) noexcept
{
	using resultType = decltype(arg1() & F2);

	field<resultType, typeOfField> retField { F1.meshRef(), resultType { 0 } };

	for (std::size_t i = 0; i < F1.size(); ++i)
		retField.ref_r()[i] = F1.ref()[i] & F2;

	return retField;
}

template<typename arg1, typename arg2, typename typeOfField>
auto dampProduct(const field<arg1, typeOfField> & F1,
		const field<arg2, typeOfField> & F2) noexcept
{
	using resultType = decltype(arg1() && arg2());

	field<resultType, typeOfField> retField { F1.meshRef(), resultType { 0 } };

	for (std::size_t i = 0; i < F1.size(); ++i)
		retField.ref_r()[i] = F1.ref()[i] && F2.ref()[i];

	return retField;
}

template<typename arg1, typename arg2, typename typeOfField>
auto dampProduct(const field<arg1, typeOfField> & F1, const arg2 & F2) noexcept
{
	using resultType = decltype(arg1() && F2);

	field<resultType, typeOfField> retField { F1.meshRef(), resultType { 0 } };

	for (std::size_t i = 0; i < F1.size(); ++i)
		retField.ref_r()[i] = F1.ref()[i] && F2;

	return retField;
}

template<typename arg1, typename arg2, typename typeOfField>
auto division(const field<arg1, typeOfField> & F1,
		const field<arg2, typeOfField> & F2) noexcept
{
	using resultType = decltype(arg1() / arg2());

	field<resultType, typeOfField> retField { F1.meshRef(), resultType { 0 } };

	for (std::size_t i = 0; i < F1.size(); ++i)
		retField.ref_r()[i] = F1.ref()[i] / F2.ref()[i];

	return retField;
}

template<typename arg1, typename arg2, typename typeOfField>
auto division(const field<arg1, typeOfField> & F1, const arg2 & F2) noexcept
{
	using resultType = decltype(arg1() / F2);

	field<resultType, typeOfField> retField { F1.meshRef(), resultType { 0 } };

	for (std::size_t i = 0; i < F1.size(); ++i)
		retField.ref_r()[i] = F1.ref()[i] / F2;

	return retField;
}

template<typename arg1, typename arg2, typename typeOfField>
auto& divisionSelf(field<arg1, typeOfField> & F1,
		const field<arg2, typeOfField> & F2) noexcept
{
	for (std::size_t i = 0; i < F1.size(); ++i)
		F1.ref_r()[i] /= F2.ref()[i];

	return F1;
}

template<typename arg1, typename arg2, typename typeOfField>
auto& divisionSelf(field<arg1, typeOfField> & F1, const arg2 & F2) noexcept
{
	for (std::size_t i = 0; i < F1.size(); ++i)
		F1.ref_r()[i] /= F2;

	return F1;
}
}  // namespace schemi

#endif /* FIELDPRODUCTS_HPP_ */
