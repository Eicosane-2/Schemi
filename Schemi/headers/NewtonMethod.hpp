/*
 * NewtonMethod.hpp
 *
 *  Created on: 2024/02/27
 *      Author: Maxim Boldyrev
 */

#ifndef NEWTONMETHOD_HPP_
#define NEWTONMETHOD_HPP_

#include <functional>

#include "scalar.hpp"

namespace schemi
{
scalar NewtonMethod(const scalar startingValue,
		const std::function<scalar(const scalar)> & f,
		const std::function<scalar(const scalar)> & dfdx) noexcept;
}  // namespace schemi

#endif /* NEWTONMETHOD_HPP_ */
