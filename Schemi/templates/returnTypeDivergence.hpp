/*
 * returnTypeDivergence.hpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#ifndef RETURNTYPEDIVERGENCE_HPP_
#define RETURNTYPEDIVERGENCE_HPP_

#include "vectorVectorDotProduct.hpp"

namespace schemi
{
template<typename T>
using returnTypeDivergence = decltype(T() & vector());
}  // namespace schemi

#endif /* RETURNTYPEDIVERGENCE_HPP_ */
