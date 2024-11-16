/*
 * zoneUnderPlanes.cpp
 *
 *  Created on: 2024/10/09
 *      Author: Maxim Boldyrev
 */

#include "zoneUnderPlanes.hpp"

schemi::zoneUnderPlanes::zoneUnderPlanes(const std::vector<vector> & ns,
		const std::vector<vector> & ps) :
		planeList(ns.size())
{
	for (std::size_t i = 0; i < planeList.size(); ++i)
		planeList[i] = plane(ns[i], ps[i]);
}

bool schemi::zoneUnderPlanes::isPointInZone(const vector & p) const noexcept
{
	bool isPointInZone(true);

	for (std::size_t i = 0; i < planeList.size(); ++i)
	{
		const auto pointUnderPlane_i = planeList[i].isPointUnderPlane(p);

		isPointInZone = isPointInZone && pointUnderPlane_i;
	}

	return isPointInZone;
}
