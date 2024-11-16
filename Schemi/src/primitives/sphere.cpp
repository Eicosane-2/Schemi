/*
 * sphere.cpp
 *
 *  Created on: 2024/10/10
 *      Author: Maxim Boldyrev
 */

#include "sphere.hpp"

#include "exception.hpp"

schemi::sphere::sphere(const vector & c, const scalar r) :
		center(c), radius(r)
{
	if (r <= 0)
		throw exception("Radius cannot be less or equal zero.",
				errors::initialisationError);
}

bool schemi::sphere::isPointInsideSphere(const vector & p) const noexcept
{
	return (p - center).mag() <= radius;
}
