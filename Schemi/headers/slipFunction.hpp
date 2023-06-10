/*
 * slipFunction.hpp
 *
 *  Created on: 2019/11/26
 *      Author: Maxim Boldyrev
 *
 *      Functions returning value for slip boundary condition.
 */

#ifndef SLIPFUNCTION_HPP_
#define SLIPFUNCTION_HPP_

#include "scalar.hpp"
#include "tensor.hpp"
#include "tensor3.hpp"
#include "vector.hpp"

namespace schemi
{
scalar slipFunction(const scalar inScalar, const vector&) noexcept;

vector slipFunction(const vector & inVector, const vector & normal) noexcept;

tensor slipFunction(const tensor & inTensor, const vector & normal) noexcept;

tensor3 slipFunction([[maybe_unused]] const tensor3 & inTensor,
		[[maybe_unused]] const vector & normal);
}  // namespace schemi

#endif /* SLIPFUNCTION_HPP_ */
