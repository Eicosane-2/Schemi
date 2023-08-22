/*
 * globalConstants.hpp
 *
 *  Created on: 2019/12/23
 *      Author: Maxim Boldyrev
 *
 *      Useful constants.
 */

#ifndef GLOBALCONSTANTS_HPP_
#define GLOBALCONSTANTS_HPP_

#include <limits>
#include <cmath>

#include "scalar.hpp"

namespace schemi
{
constexpr scalar onethirds { 1. / 3. },

twothirds { 2. / 3. },

NAvogardro { 6.02214076E23 },

R_SI { 8.31446261815324 },

PlanckConstant_SI { 6.626E-34 },

/*std::numeric_limits<scalar>::min() in denominator can lead to infinity.*/
stabilizator { 5 * std::numeric_limits<scalar>::epsilon() },

veryBig { 1E-3 * std::numeric_limits<scalar>::max() },

#if ( defined(__GNUG__) ) && ( !defined(__ICC) )
		Pi_number { std::acos(-1.0) }, e_number { std::exp(1.0) },
#else
		Pi_number { 3.1415926535897932384626433832795 },
		e_number { 2.71828182845904523536 },
#endif

		zeroLevel { stabilizator },

		zeroMix { 0.01 },

		kappaPPM { onethirds };

constexpr std::size_t componentPlaceholder { static_cast<std::size_t>(-1) },
		lengthOfNumber { 6 };
constexpr int ioPrecision { 20 };
}  // namespace schemi

#endif /* GLOBALCONSTANTS_HPP_ */
