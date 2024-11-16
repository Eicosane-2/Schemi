/*
 * zone.hpp
 *
 *  Created on: 2024/10/09
 *      Author: Maxim Boldyrev
 */

#ifndef ZONE_HPP_
#define ZONE_HPP_

#include <vector>
#include <memory>

#include "vector.hpp"

namespace schemi
{
class zone
{
public:
	virtual ~zone() noexcept =0;

	virtual bool isPointInZone(const vector & p) const noexcept =0;

	static std::vector<std::unique_ptr<zone>> zonesArray();
};
}  // namespace schemi

#endif /* ZONE_HPP_ */
