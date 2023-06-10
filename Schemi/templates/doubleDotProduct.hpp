/*
 * doubleDotProduct.hpp
 *
 *  Created on: 2019/11/13
 *      Author: Maxim Boldyrev
 *
 *      Function for double dot product calculation.
 */

#ifndef DOUBLEDOTPRODUCT_HPP_
#define DOUBLEDOTPRODUCT_HPP_

namespace schemi
{
template<typename T>
scalar operator&&(const T & inTensor1, const T & inTensor2) noexcept
{
	scalar result { 0 };

	for (std::size_t i = 0; i < T::vsize; ++i)
		result += inTensor1.v()[i] * inTensor2.v()[i];

	return result;
}
}  // namespace schemi

#endif /* DOUBLEDOTPRODUCT_HPP_ */
