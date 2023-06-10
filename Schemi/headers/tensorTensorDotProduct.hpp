/*
 * tensorTensorDotProduct.hpp
 *
 *  Created on: 2021/05/29
 *      Author: Maxim Boldyrev
 *
 *      Function for 2nd rank tensor-tensor dot product.
 */

#ifndef TENSORTENSORDOTPRODUCT_HPP_
#define TENSORTENSORDOTPRODUCT_HPP_

#include "tensor.hpp"

namespace schemi
{
tensor operator&(const tensor & inTensor1, const tensor & inTensor2) noexcept;
}

#endif /* TENSORTENSORDOTPRODUCT_HPP_ */
