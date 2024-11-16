/*
 * zoneInsideSphere.hpp
 *
 *  Created on: 2024/10/10
 *      Author: Maxim Boldyrev
 */

#ifndef ZONEINSIDESPHERE_HPP_
#define ZONEINSIDESPHERE_HPP_

#include "zone.hpp"

#include "sphere.hpp"

namespace schemi
{
class zoneInsideSphere: public zone, private sphere
{
public:
	zoneInsideSphere(const vector & c, const scalar r);

	bool isPointInZone(const vector & p) const noexcept override;
};
}

#endif /* ZONEINSIDESPHERE_HPP_ */
