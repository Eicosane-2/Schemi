/*
 * noslipFunction.hpp
 *
 *  Created on: 2026/09/08
 *      Author: Maxim Boldyrev
 *
 *      Functions returning value for no-slip boundary condition.
 */

#ifndef NOSLIPFUNCTION_HPP_
#define NOSLIPFUNCTION_HPP_

#include "scalar.hpp"
#include "tensor.hpp"
#include "tensor3.hpp"
#include "vector.hpp"

namespace schemi
{
scalar noslipFunction(const scalar inScalar, const vector&) noexcept;

vector noslipFunction(const vector&, const vector&) noexcept;

tensor noslipFunction(const tensor&, const vector&) noexcept;

tensor3 noslipFunction(const tensor3&, const vector&);
}  // namespace schemi

#endif /* NOSLIPFUNCTION_HPP_ */
