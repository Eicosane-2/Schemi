/*
 * zoneUnderPlane.cpp
 *
 *  Created on: 2024/10/09
 *      Author: Maxim Boldyrev
 */

#include "zoneUnderPlane.hpp"

schemi::zoneUnderPlane::zoneUnderPlane(const vector & n, const vector & p) :
		plane(n, p)
{
}

bool schemi::zoneUnderPlane::isPointInZone(const vector & p) const noexcept
{
	return isPointUnderPlane(p);
}
