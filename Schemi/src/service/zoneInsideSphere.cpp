/*
 * zoneInsideSphere.cpp
 *
 *  Created on: 2024/10/10
 *      Author: Maxim Boldyrev
 */

#include "zoneInsideSphere.hpp"

schemi::zoneInsideSphere::zoneInsideSphere(const vector & c, const scalar r) :
		sphere(c, r)
{
}

bool schemi::zoneInsideSphere::isPointInZone(const vector & p) const noexcept
{
	return isPointInsideSphere(p);
}
