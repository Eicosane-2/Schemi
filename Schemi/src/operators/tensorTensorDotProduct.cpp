/*
 * tensorTensorDotProduct.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensorTensorDotProduct.hpp"

schemi::tensor schemi::operator&(const tensor & inTensor1,
		const tensor & inTensor2) noexcept
{
	return tensor { inTensor1.v()[0] * inTensor2.v()[0]
			+ inTensor1.v()[1] * inTensor2.v()[3]
			+ inTensor1.v()[2] * inTensor2.v()[6], inTensor1.v()[0]
			* inTensor2.v()[1] + inTensor1.v()[1] * inTensor2.v()[4]
			+ inTensor1.v()[2] * inTensor2.v()[7], inTensor1.v()[0]
			* inTensor2.v()[2] + inTensor1.v()[1] * inTensor2.v()[5]
			+ inTensor1.v()[2] * inTensor2.v()[8], inTensor1.v()[3]
			* inTensor2.v()[0] + inTensor1.v()[4] * inTensor2.v()[3]
			+ inTensor1.v()[5] * inTensor2.v()[6], inTensor1.v()[3]
			* inTensor2.v()[1] + inTensor1.v()[4] * inTensor2.v()[4]
			+ inTensor1.v()[5] * inTensor2.v()[7], inTensor1.v()[3]
			* inTensor2.v()[2] + inTensor1.v()[4] * inTensor2.v()[5]
			+ inTensor1.v()[5] * inTensor2.v()[8], inTensor1.v()[6]
			* inTensor2.v()[0] + inTensor1.v()[7] * inTensor2.v()[3]
			+ inTensor1.v()[8] * inTensor2.v()[6], inTensor1.v()[6]
			* inTensor2.v()[1] + inTensor1.v()[7] * inTensor2.v()[4]
			+ inTensor1.v()[8] * inTensor2.v()[7], inTensor1.v()[6]
			* inTensor2.v()[2] + inTensor1.v()[7] * inTensor2.v()[5]
			+ inTensor1.v()[8] * inTensor2.v()[8] };
}
