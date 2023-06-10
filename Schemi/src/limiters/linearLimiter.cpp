/*
 * linearLimiter.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "linearLimiter.hpp"

schemi::vector schemi::linearLimiter::calculate(const vector&,
		const vector & gradientC) const noexcept
{
	return gradientC;
}

schemi::tensor schemi::linearLimiter::calculate(const tensor&,
		const tensor & gradientC) const noexcept
{
	return gradientC;
}

schemi::tensor3 schemi::linearLimiter::calculate(const tensor3&,
		const tensor3 & gradientC) const noexcept
{
	return gradientC;
}

schemi::vector schemi::linearLimiter::calculateNoRightLimit(const vector&,
		const vector & gradientC) const noexcept
{
	return gradientC;
}

schemi::tensor schemi::linearLimiter::calculateNoRightLimit(const tensor&,
		const tensor & gradientC) const noexcept
{
	return gradientC;
}

schemi::tensor3 schemi::linearLimiter::calculateNoRightLimit(const tensor3&,
		const tensor3 & gradientC) const noexcept
{
	return gradientC;
}

