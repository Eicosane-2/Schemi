/*
 * vectorVectorDotProduct.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vectorVectorDotProduct.hpp"

#include <cstddef>

#include "scalar.hpp"
#include "vector.hpp"

schemi::scalar schemi::operator&(const schemi::vector & inVector1,
		const schemi::vector & inVector2) noexcept
{
	scalar result = 0;

	for (std::size_t i = 0; i < vector::vsize; ++i)
		result += inVector1.v()[i] * inVector2.v()[i];

	return result;
}
