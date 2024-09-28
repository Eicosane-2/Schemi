/*
 * secantMethod.cpp
 *
 *  Created on: 2024/02/18
 *      Author: Maxim Boldyrev
 */

#include "secantMethod.hpp"

#include "exception.hpp"
#include "globalConstants.hpp"

schemi::scalar schemi::secantMethod(const scalar startingValue,
		const scalar guess, const std::function<scalar(const scalar)> & f)
{
	constexpr static std::size_t maxIter { 100 };

	std::size_t it { 0 };

	scalar x1 { guess }, x2 { startingValue };

	if (std::abs(x2 - x1) <= convergenceToleranceGlobal)
		throw exception("Initial values too close.",
				errors::initialisationError);

	while ((std::abs(x2 - x1) > convergenceToleranceGlobal) && (it <= maxIter))
	{
		const auto x3 = x2 - f(x2) * (x2 - x1) / (f(x2) - f(x1));

		x1 = x2;
		x2 = x3;

		it++;
	}

	return x2;
}
