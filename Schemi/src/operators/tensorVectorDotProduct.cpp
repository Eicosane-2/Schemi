/*
 * tensorVectorDotProduct.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensorVectorDotProduct.hpp"

schemi::vector schemi::operator&(const tensor & inTensor,
		const vector & inVector) noexcept
{
	return vector { inTensor.v()[0] * inVector.v()[0]
			+ inTensor.v()[1] * inVector.v()[1]
			+ inTensor.v()[2] * inVector.v()[2], inTensor.v()[3]
			* inVector.v()[0] + inTensor.v()[4] * inVector.v()[1]
			+ inTensor.v()[5] * inVector.v()[2], inTensor.v()[6]
			* inVector.v()[0] + inTensor.v()[7] * inVector.v()[1]
			+ inTensor.v()[8] * inVector.v()[2] };
}

schemi::vector schemi::operator&(const vector & inVector,
		const tensor & inTensor) noexcept
{
	tensor bufTensor(inTensor);
	bufTensor.transpose();

	return bufTensor & inVector;
}

schemi::tensor schemi::operator&(const tensor3 & inTensor,
		const vector & inVector) noexcept
{
	return tensor { inTensor.v()[0] * inVector.v()[0]
			+ inTensor.v()[1] * inVector.v()[1]
			+ inTensor.v()[2] * inVector.v()[2], inTensor.v()[3]
			* inVector.v()[0] + inTensor.v()[4] * inVector.v()[1]
			+ inTensor.v()[5] * inVector.v()[2], inTensor.v()[6]
			* inVector.v()[0] + inTensor.v()[7] * inVector.v()[1]
			+ inTensor.v()[8] * inVector.v()[2],

	inTensor.v()[9 + 0] * inVector.v()[0]
			+ inTensor.v()[9 + 1] * inVector.v()[1]
			+ inTensor.v()[9 + 2] * inVector.v()[2], inTensor.v()[9 + 3]
			* inVector.v()[0] + inTensor.v()[9 + 4] * inVector.v()[1]
			+ inTensor.v()[9 + 5] * inVector.v()[2], inTensor.v()[9 + 6]
			* inVector.v()[0] + inTensor.v()[9 + 7] * inVector.v()[1]
			+ inTensor.v()[9 + 8] * inVector.v()[2],

	inTensor.v()[18 + 0] * inVector.v()[0]
			+ inTensor.v()[18 + 1] * inVector.v()[1]
			+ inTensor.v()[18 + 2] * inVector.v()[2], inTensor.v()[18 + 3]
			* inVector.v()[0] + inTensor.v()[18 + 4] * inVector.v()[1]
			+ inTensor.v()[18 + 5] * inVector.v()[2], inTensor.v()[18 + 6]
			* inVector.v()[0] + inTensor.v()[18 + 7] * inVector.v()[1]
			+ inTensor.v()[18 + 8] * inVector.v()[2] };
}

schemi::tensor schemi::operator&(const vector & inVector,
		const tensor3 & inTensor) noexcept
{
	return tensor { inTensor.v()[0] * inVector.v()[0]
			+ inTensor.v()[3] * inVector.v()[1]
			+ inTensor.v()[6] * inVector.v()[2], inTensor.v()[1]
			* inVector.v()[0] + inTensor.v()[4] * inVector.v()[1]
			+ inTensor.v()[7] * inVector.v()[2], inTensor.v()[2]
			* inVector.v()[0] + inTensor.v()[5] * inVector.v()[1]
			+ inTensor.v()[8] * inVector.v()[2],

	inTensor.v()[9 + 0] * inVector.v()[0]
			+ inTensor.v()[9 + 3] * inVector.v()[1]
			+ inTensor.v()[9 + 6] * inVector.v()[2], inTensor.v()[9 + 1]
			* inVector.v()[0] + inTensor.v()[9 + 4] * inVector.v()[1]
			+ inTensor.v()[9 + 7] * inVector.v()[2], inTensor.v()[9 + 2]
			* inVector.v()[0] + inTensor.v()[9 + 5] * inVector.v()[1]
			+ inTensor.v()[9 + 8] * inVector.v()[2],

	inTensor.v()[18 + 0] * inVector.v()[0]
			+ inTensor.v()[18 + 3] * inVector.v()[1]
			+ inTensor.v()[18 + 6] * inVector.v()[2], inTensor.v()[18 + 1]
			* inVector.v()[0] + inTensor.v()[18 + 4] * inVector.v()[1]
			+ inTensor.v()[18 + 7] * inVector.v()[2], inTensor.v()[18 + 2]
			* inVector.v()[0] + inTensor.v()[18 + 5] * inVector.v()[1]
			+ inTensor.v()[18 + 8] * inVector.v()[2] };
}
