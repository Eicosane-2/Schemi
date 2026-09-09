/*
 * noslipFunction.cpp
 *
 *  Created on: 2026/09/08
 *      Author: Maxim Boldyrev
 */

#include "noslipFunction.hpp"

#include "exception.hpp"
#include "tensor.hpp"
#include "vector.hpp"
#include "doubleDotProduct.hpp"

schemi::scalar schemi::noslipFunction(const scalar inScalar,
		const vector&) noexcept
{
	return inScalar;
}

schemi::vector schemi::noslipFunction(const vector&, const vector&) noexcept
{
	return vector(0);
}

schemi::tensor schemi::noslipFunction(const tensor&, const vector&) noexcept
{
	return tensor(0);
}

schemi::tensor3 schemi::noslipFunction(const tensor3&, const vector&)
{
	throw exception(
			"<<noslipFunction>> is not implemented for a third rank tensor (tensor3).",
			errors::boundaryConditionError);
}
