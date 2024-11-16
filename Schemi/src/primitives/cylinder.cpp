/*
 * cylinder.cpp
 *
 *  Created on: 2024/10/10
 *      Author: Maxim Boldyrev
 */

#include "cylinder.hpp"

#include "exception.hpp"

schemi::cylinder::cylinder(const vector & n1, const vector & c1,
		const vector & c2, const scalar r) :
		radius(r)
{
	const auto n2 = n1 * (-1);

	planes = { plane(n1 / n1.mag(), c1), plane(n2 / n2.mag(), c2) };

	centerLine = line(n1 / n1.mag(), c2);

	if (r <= 0)
		throw exception("Radius cannot be less or equal zero.",
				errors::initialisationError);
}

bool schemi::cylinder::isPointInsideCylinder(const vector & p) const noexcept
{
	bool isPointInsideCylinder(true);

	for (std::size_t i = 0; i < planes.size(); ++i)
	{
		const auto pointUnderPlane_i = planes[i].isPointUnderPlane(p);

		isPointInsideCylinder = isPointInsideCylinder && pointUnderPlane_i;
	}

	isPointInsideCylinder = isPointInsideCylinder
			&& (centerLine.distanceFromPoint(p) <= radius);

	return isPointInsideCylinder;
}
