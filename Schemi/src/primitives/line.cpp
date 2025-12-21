/*
 * line.cpp
 *
 *  Created on: 2024/10/10
 *      Author: Maxim Boldyrev
 */

#include "line.hpp"

#include "vector.hpp"

schemi::line::line(const vector & t, const vector & p) noexcept :
		tangent(t / t.mag()), point(p)
{
}

schemi::scalar schemi::line::distanceFromPoint(const vector & p) const noexcept
{
	return ((point - p) - ((point - p) & tangent) * tangent).mag();
}
