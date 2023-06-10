/*
 * dyadicProduct.hpp
 *
 *  Created on: 2019/11/12
 *      Author: Maxim Boldyrev
 *
 *      Functions for vectors' dyadic products.
 */

#ifndef DYADICPRODUCT_HPP_
#define DYADICPRODUCT_HPP_

#include "tensor.hpp"
#include "tensor3.hpp"
#include "vector.hpp"

namespace schemi
{
tensor operator*(const vector & vector1, const vector & vector2) noexcept;

tensor3 operator*(const tensor & inTensor, const vector & inVector) noexcept;

tensor3 operator*(const vector & inVector, const tensor & inTensor) noexcept;
}  // namespace schemi

#endif /* DYADICPRODUCT_HPP_ */
