/*
 * zoneInsideCylinder.cpp
 *
 *  Created on: 2024/10/12
 *      Author: Maxim Boldyrev
 */

#include "zoneInsideCylinder.hpp"

schemi::zoneInsideCylinder::zoneInsideCylinder(const vector & n1,
		const vector & c1, const vector & c2, const scalar r) :
		cylinder(n1, c1, c2, r)
{
}

bool schemi::zoneInsideCylinder::isPointInZone(const vector & p) const noexcept
{
	return isPointInsideCylinder(p);
}
