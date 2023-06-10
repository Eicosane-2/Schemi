/*
 * starFields.hpp
 *
 *  Created on: 2019/12/23
 *      Author: Maxim Boldyrev
 *
 *      Structure storing some fields' values in star regions (on surfaces of mesh).
 */

#ifndef STARFIELDS_HPP_
#define STARFIELDS_HPP_

#include <vector>

#include "surfaceField.hpp"

namespace schemi
{
struct starFields
{
	std::vector<surfaceField<scalar>> c;
	surfaceField<scalar> rho;
	surfaceField<vector> v;
	surfaceField<scalar> p;
	surfaceField<vector> a;
	surfaceField<scalar> b;

	starFields(const mesh & meshRef, const std::size_t numberOfcomponents);
};
}  // namespace schemi

#endif /* STARFIELDS_HPP_ */
