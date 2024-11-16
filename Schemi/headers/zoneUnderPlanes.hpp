/*
 * zoneUnderPlanes.hpp
 *
 *  Created on: 2024/10/09
 *      Author: Maxim Boldyrev
 */

#ifndef ZONEUNDERPLANES_HPP_
#define ZONEUNDERPLANES_HPP_

#include "zone.hpp"

#include <vector>

#include "plane.hpp"

namespace schemi
{
class zoneUnderPlanes: public zone
{
	std::vector<plane> planeList;
public:
	zoneUnderPlanes(const std::vector<vector> & ns,
			const std::vector<vector> & ps);

	bool isPointInZone(const vector & p) const noexcept override;
};
}  // namespace schemi

#endif /* ZONEUNDERPLANES_HPP_ */
