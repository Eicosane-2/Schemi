/*
 * tensorVectorDotProduct.hpp
 *
 *  Created on: 2019/11/13
 *      Author: Maxim Boldyrev
 *
 *      Functions for tensor-vector dot product.
 */

#ifndef TENSORVECTORDOTPRODUCT_HPP_
#define TENSORVECTORDOTPRODUCT_HPP_

#include "tensor.hpp"
#include "tensor3.hpp"
#include "vector.hpp"

namespace schemi
{
vector operator&(const tensor & inTensor, const vector & inVector) noexcept;

vector operator&(const vector & inVector, const tensor & inTensor) noexcept;

tensor operator&(const tensor3 & inTensor, const vector & inVector) noexcept;

tensor operator&(const vector & inVector, const tensor3 & inTensor) noexcept;
}  // namespace schemi

#endif /* TENSORVECTORDOTPRODUCT_HPP_ */
