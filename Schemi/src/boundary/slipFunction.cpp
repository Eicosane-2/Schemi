/*
 * slipFunction.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "slipFunction.hpp"

#include "dyadicProduct.hpp"
#include "exception.hpp"
#include "tensorTensorDotProduct.hpp"
#include "tensorVectorDotProduct.hpp"
#include "vectorVectorDotProduct.hpp"
#include "doubleDotProduct.hpp"

schemi::scalar schemi::slipFunction(const scalar inScalar,
		const vector&) noexcept
{
	return inScalar;
}

schemi::vector schemi::slipFunction(const vector & inVector,
		const vector & normal) noexcept
{
	return inVector - (inVector & normal) * normal;
}

schemi::tensor schemi::slipFunction(const tensor & inTensor,
		const vector & normal) noexcept
{
	const tensor normalTensor(normal * normal);

	return inTensor
			- ((inTensor & normalTensor) + (normalTensor & inTensor)
					- 2 * (inTensor && normalTensor) * normalTensor);
}

schemi::tensor3 schemi::slipFunction([[maybe_unused]] const tensor3 & inTensor,
		[[maybe_unused]] const vector & normal)
{
	throw exception(
			"<<slipFunction>> is not implemented for a third rank tensor (tensor3).",
			errors::boundaryConditionError);

	return tensor3(0);
}
