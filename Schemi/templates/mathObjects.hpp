/*
 * mathObjects.hpp
 *
 *  Created on: 2025/11/22
 *      Author: Maxim Boldyrev
 */

#ifndef MATHOBJECTS_HPP_
#define MATHOBJECTS_HPP_

#include <concepts>

#include "scalar.hpp"
#include "tensor.hpp"
#include "tensor3.hpp"
#include "vector.hpp"

namespace schemi
{
template<typename T>
concept mathObjects = std::same_as<T, int> || std::same_as<T, scalar> ||
std::same_as<T, vector> || std::same_as<T, tensor> ||
std::same_as<T, tensor3>;

template<typename T>
concept tensors = std::same_as<T, tensor> || std::same_as<T, tensor3>;
} // namespace schemi

#endif /* MATHOBJECTS_HPP_ */
