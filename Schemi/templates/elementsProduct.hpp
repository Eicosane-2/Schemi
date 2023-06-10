/*
 * elementsProduct.hpp
 *
 *  Created on: 2022/12/04
 *      Author: Maxim Boldyrev
 */

#ifndef ELEMENTSPRODUCT_HPP_
#define ELEMENTSPRODUCT_HPP_

namespace schemi
{
template<typename T>
inline T elementsProduct(const T & a, const T & b) noexcept
{
	T ret;

	for (std::size_t i = 0; i < T::vsize; ++i)
		ret.v_r()[i] = a.v()[i] * b.v()[i];

	return ret;
}

template<typename T>
inline T elementsDivision(const T & a, const T & b) noexcept
{
	T ret;

	for (std::size_t i = 0; i < T::vsize; ++i)
		ret.v_r()[i] = a.v()[i] / (b.v()[i] + stabilizator);

	return ret;
}
}  // namespace schemi

#endif /* ELEMENTSPRODUCT_HPP_ */
