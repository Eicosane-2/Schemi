/*
 * dyadicProduct.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "dyadicProduct.hpp"

schemi::tensor schemi::operator*(const vector & vector1,
		const vector & vector2) noexcept
{
	return tensor { vector1.v()[0] * vector2.v()[0], vector1.v()[0]
			* vector2.v()[1], vector1.v()[0] * vector2.v()[2], vector1.v()[1]
			* vector2.v()[0], vector1.v()[1] * vector2.v()[1], vector1.v()[1]
			* vector2.v()[2], vector1.v()[2] * vector2.v()[0], vector1.v()[2]
			* vector2.v()[1], vector1.v()[2] * vector2.v()[2] };
}

schemi::tensor3 schemi::operator*(const tensor & inTensor,
		const vector & inVector) noexcept
{
	return tensor3 {

	inTensor.v()[0] * inVector.v()[0], inTensor.v()[0] * inVector.v()[1],
			inTensor.v()[0] * inVector.v()[2], inTensor.v()[1]
					* inVector.v()[0], inTensor.v()[1] * inVector.v()[1],
			inTensor.v()[1] * inVector.v()[2], inTensor.v()[2]
					* inVector.v()[0], inTensor.v()[2] * inVector.v()[1],
			inTensor.v()[2] * inVector.v()[2],

			inTensor.v()[3] * inVector.v()[0], inTensor.v()[3]
					* inVector.v()[1], inTensor.v()[3] * inVector.v()[2],
			inTensor.v()[4] * inVector.v()[0], inTensor.v()[4]
					* inVector.v()[1], inTensor.v()[4] * inVector.v()[2],
			inTensor.v()[5] * inVector.v()[0], inTensor.v()[5]
					* inVector.v()[1], inTensor.v()[5] * inVector.v()[2],

			inTensor.v()[6] * inVector.v()[0], inTensor.v()[6]
					* inVector.v()[1], inTensor.v()[6] * inVector.v()[2],
			inTensor.v()[7] * inVector.v()[0], inTensor.v()[7]
					* inVector.v()[1], inTensor.v()[7] * inVector.v()[2],
			inTensor.v()[8] * inVector.v()[0], inTensor.v()[8]
					* inVector.v()[1], inTensor.v()[8] * inVector.v()[2] };
}

schemi::tensor3 schemi::operator*(const vector & inVector,
		const tensor & inTensor) noexcept
{
	return tensor3 {

	inVector.v()[0] * inTensor.v()[0], inVector.v()[0] * inTensor.v()[1],
			inVector.v()[0] * inTensor.v()[2], inVector.v()[0]
					* inTensor.v()[3], inVector.v()[0] * inTensor.v()[4],
			inVector.v()[0] * inTensor.v()[5], inVector.v()[0]
					* inTensor.v()[6], inVector.v()[0] * inTensor.v()[7],
			inVector.v()[0] * inTensor.v()[8],

			inVector.v()[1] * inTensor.v()[0], inVector.v()[1]
					* inTensor.v()[1], inVector.v()[1] * inTensor.v()[2],
			inVector.v()[1] * inTensor.v()[3], inVector.v()[1]
					* inTensor.v()[4], inVector.v()[1] * inTensor.v()[5],
			inVector.v()[1] * inTensor.v()[6], inVector.v()[1]
					* inTensor.v()[7], inVector.v()[1] * inTensor.v()[8],

			inVector.v()[2] * inTensor.v()[0], inVector.v()[2]
					* inTensor.v()[1], inVector.v()[2] * inTensor.v()[2],
			inVector.v()[2] * inTensor.v()[3], inVector.v()[2]
					* inTensor.v()[4], inVector.v()[2] * inTensor.v()[5],
			inVector.v()[2] * inTensor.v()[6], inVector.v()[2]
					* inTensor.v()[7], inVector.v()[2] * inTensor.v()[8] };
}
