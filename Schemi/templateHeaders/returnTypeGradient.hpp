/*
 * returnTypeGradient.hpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#ifndef RETURNTYPEGRADIENT_HPP_
#define RETURNTYPEGRADIENT_HPP_

namespace schemi
{
template<typename T>
using returnTypeGradient = decltype(T() * vector());
}  // namespace schemi

#endif /* RETURNTYPEGRADIENT_HPP_ */
