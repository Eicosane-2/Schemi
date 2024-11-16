/*
 * plane.hpp
 *
 *  Created on: 2024/10/09
 *      Author: Maxim Boldyrev
 */

#ifndef PLANE_HPP_
#define PLANE_HPP_

#include "vector.hpp"

namespace schemi
{
class plane
{
	vector normale { 0, 0, 0 };
	vector point { 0, 0, 0 };
public:
	plane() noexcept = default;

	plane(const vector & n, const vector & p) noexcept;

	bool isPointUnderPlane(const vector & p) const noexcept;
};
}

#endif /* PLANE_HPP_ */
