/*
 * zoneInsideCylinder.hpp
 *
 *  Created on: 2024/10/12
 *      Author: Maxim Boldyrev
 */

#ifndef ZONEINSIDECYLINDER_HPP_
#define ZONEINSIDECYLINDER_HPP_

#include "zone.hpp"

#include "cylinder.hpp"

namespace schemi
{
class zoneInsideCylinder: public zone, private cylinder
{
public:
	zoneInsideCylinder(const vector & n1, const vector & c1, const vector & c2,
			const scalar r);

	bool isPointInZone(const vector & p) const noexcept override;
};
}

#endif /* ZONEINSIDECYLINDER_HPP_ */
