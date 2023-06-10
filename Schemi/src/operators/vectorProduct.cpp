/*
 * vectorProduct.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vectorProduct.hpp"

schemi::vector schemi::vectorProduct(const schemi::vector & a,
		const schemi::vector & b) noexcept
{
	return vector { a.v()[1] * b.v()[2] - a.v()[2] * b.v()[1], -(a.v()[0]
			* b.v()[2] - a.v()[2] * b.v()[0]), a.v()[0] * b.v()[1]
			- a.v()[1] * b.v()[0] };
}
