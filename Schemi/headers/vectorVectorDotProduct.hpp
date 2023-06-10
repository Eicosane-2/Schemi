/*
 * vectorVectorDotProduct.hpp
 *
 *  Created on: 2019/11/13
 *      Author: Maxim Boldyrev
 *
 *      Function for vector-vector dot product.
 */

#ifndef VECTORVECTORDOTPRODUCT_HPP_
#define VECTORVECTORDOTPRODUCT_HPP_

#include "scalar.hpp"
#include "vector.hpp"

namespace schemi
{
scalar operator&(const vector & inVector1, const vector & inVector2) noexcept;
}  // namespace schemi

#endif /* VECTORVECTORDOTPRODUCT_HPP_ */
