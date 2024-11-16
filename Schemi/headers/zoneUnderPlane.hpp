/*
 * zoneUnderPlane.hpp
 *
 *  Created on: 2024/10/09
 *      Author: Maxim Boldyrev
 */

#ifndef ZONEUNDERPLANE_HPP_
#define ZONEUNDERPLANE_HPP_

#include "zone.hpp"

#include "plane.hpp"

namespace schemi
{
class zoneUnderPlane: public zone, private plane
{
public:
	zoneUnderPlane(const vector & n, const vector & p);

	bool isPointInZone(const vector & p) const noexcept override;
};
}

#endif /* ZONEUNDERPLANE_HPP_ */
