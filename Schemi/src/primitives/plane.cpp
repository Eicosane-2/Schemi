/*
 * plane.cpp
 *
 *  Created on: 2024/10/09
 *      Author: Maxim Boldyrev
 */

#include "plane.hpp"

#include "vector.hpp"

schemi::plane::plane(const vector & n, const vector & p) noexcept :
		normale(n / n.mag()), point(p)
{
}

bool schemi::plane::isPointUnderPlane(const vector & p) const noexcept
{
	return ((point - p) & normale) >= 0;
}
