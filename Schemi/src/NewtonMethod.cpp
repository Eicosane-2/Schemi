/*
 * NewtonMethod.cpp
 *
 *  Created on: 2024/02/27
 *      Author: Maxim Boldyrev
 */

#include "NewtonMethod.hpp"

#include "exception.hpp"
#include "globalConstants.hpp"

schemi::scalar schemi::NewtonMethod(const scalar startingValue,
		const std::function<scalar(const scalar)> & f,
		const std::function<scalar(const scalar)> & dfdx) noexcept
{
	constexpr static std::size_t maxIter { 100 };

	std::size_t it { 1 };

	//First iteration.
	scalar x1 { startingValue };

	scalar x2 = x1 - f(x1) / dfdx(x1);

	if (std::abs(x2 - x1) <= convergenceToleranceGlobal)
		return x2;

	//Next iterations.
	while ((std::abs(x2 - x1) > convergenceToleranceGlobal) && (it <= maxIter))
	{
		x1 = x2;

		x2 = x1 - f(x1) / dfdx(x1);

		it++;
	}

	return x2;
}
