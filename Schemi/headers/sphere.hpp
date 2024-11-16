/*
 * sphere.hpp
 *
 *  Created on: 2024/10/10
 *      Author: Maxim Boldyrev
 */

#ifndef SPHERE_HPP_
#define SPHERE_HPP_

#include "vector.hpp"

namespace schemi
{
class sphere
{
	vector center { 0, 0, 0 };
	scalar radius { 0 };
public:
	sphere() noexcept = default;

	sphere(const vector & c, const scalar r);

	bool isPointInsideSphere(const vector & p) const noexcept;
};
}

#endif /* SPHERE_HPP_ */
