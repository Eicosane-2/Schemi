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
	return vector { std::get<1>(a()) * std::get<2>(b())
			- std::get<2>(a()) * std::get<1>(b()), -(std::get<0>(a())
			* std::get<2>(b()) - std::get<2>(a()) * std::get<0>(b())), std::get<
			0>(a()) * std::get<1>(b()) - std::get<1>(a()) * std::get<0>(b()) };
}
