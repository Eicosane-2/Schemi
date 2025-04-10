/*
 * linearLimiter.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "linearLimiter.hpp"

schemi::vector schemi::linearLimiter::calculate(const vector&,
		const vector & gradient) const noexcept
{
	return gradient;
}

schemi::tensor schemi::linearLimiter::calculate(const tensor&,
		const tensor & gradient) const noexcept
{
	return gradient;
}

schemi::tensor3 schemi::linearLimiter::calculate(const tensor3&,
		const tensor3 & gradient) const noexcept
{
	return gradient;
}

schemi::vector schemi::linearLimiter::calculateNoRSLimit(const vector&,
		const vector & gradient) const noexcept
{
	return gradient;
}

schemi::tensor schemi::linearLimiter::calculateNoRSLimit(const tensor&,
		const tensor & gradient) const noexcept
{
	return gradient;
}

schemi::tensor3 schemi::linearLimiter::calculateNoRSLimit(const tensor3&,
		const tensor3 & gradient) const noexcept
{
	return gradient;
}

