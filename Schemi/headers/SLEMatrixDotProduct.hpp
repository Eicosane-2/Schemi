/*
 * SLEMatrixDotProduct.hpp
 *
 *  Created on: 2021/01/23
 *      Author: Maxim Boldyrev
 *
 *      Function for matrix-vector dot product.
 */

#ifndef SLEMATRIXDOTPRODUCT_HPP_
#define SLEMATRIXDOTPRODUCT_HPP_

#include <valarray>

#include "scalar.hpp"
#include "SLEMatrix.hpp"

namespace schemi
{
std::valarray<scalar> operator&(const SLEMatrix::SLEMatrixStorage & M,
		const std::valarray<scalar> & v) noexcept;
}  // namespace schemi

#endif /* SLEMATRIXDOTPRODUCT_HPP_ */
