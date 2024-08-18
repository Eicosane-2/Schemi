/*
 * zeroLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "zeroLimiter.hpp"

schemi::vector schemi::zeroLimiter::calculate(const vector&,
		const vector&) const noexcept
{
	return vector { 0 };
}

schemi::tensor schemi::zeroLimiter::calculate(const tensor&,
		const tensor&) const noexcept
{
	return tensor { 0 };
}

schemi::tensor3 schemi::zeroLimiter::calculate(const tensor3&,
		const tensor3&) const noexcept
{
	return tensor3 { 0 };
}

schemi::vector schemi::zeroLimiter::calculateNoRSLimit(const vector&,
		const vector&) const noexcept
{
	return vector { 0 };
}

schemi::tensor schemi::zeroLimiter::calculateNoRSLimit(const tensor&,
		const tensor&) const noexcept
{
	return tensor { 0 };
}

schemi::tensor3 schemi::zeroLimiter::calculateNoRSLimit(const tensor3&,
		const tensor3&) const noexcept
{
	return tensor3 { 0 };
}
