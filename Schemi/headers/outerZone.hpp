/*
 * outerZone.hpp
 *
 *  Created on: 2024/10/09
 *      Author: Maxim Boldyrev
 */

#ifndef OUTERZONE_HPP_
#define OUTERZONE_HPP_

#include "zone.hpp"

namespace schemi
{
class outerZone: public zone
{
	bool isPointInZone(const vector & p) const noexcept override;
};
}  // namespace schemi

#endif /* OUTERZONE_HPP_ */
