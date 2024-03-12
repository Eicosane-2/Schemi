/*
 * secantMethod.hpp
 *
 *  Created on: 2024/02/18
 *      Author: Maxim Boldyrev
 */

#ifndef SECANTMETHOD_HPP_
#define SECANTMETHOD_HPP_

#include <functional>

#include "scalar.hpp"

namespace schemi
{
scalar secantMethod(const scalar startingValue, const scalar guess,
		const std::function<scalar(const scalar)> & f);
}  // namespace schemi

#endif /* SECANTMETHOD_HPP_ */
