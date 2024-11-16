/*
 * cylinder.hpp
 *
 *  Created on: 2024/10/10
 *      Author: Maxim Boldyrev
 */

#ifndef CYLINDER_HPP_
#define CYLINDER_HPP_

#include "line.hpp"
#include "plane.hpp"
#include "vector.hpp"

namespace schemi
{
class cylinder: public line
{
	std::array<plane, 2> planes { plane(), plane() };
	scalar radius { 0 };
	line centerLine { vector(0), vector(0) };
public:
	cylinder() noexcept = default;

	cylinder(const vector & n1, const vector & c1, const vector & c2,
			const scalar r);

	bool isPointInsideCylinder(const vector & p) const noexcept;
};
}

#endif /* CYLINDER_HPP_ */
