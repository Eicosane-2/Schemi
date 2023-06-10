/*
 * surfaceField.hpp
 *
 *  Created on: 2019/11/23
 *      Author: Maxim Boldyrev
 *
 *      Alias template for surface field.
 */

#ifndef SURFACEFIELD_HPP_
#define SURFACEFIELD_HPP_

#include "quadraticSurface.hpp"
#include "field.hpp"

namespace schemi
{
template<typename typeOfValue> using surfaceField = field<typeOfValue, quadraticSurface>;
}  // namespace schemi

#endif /* SURFACEFIELD_HPP_ */
