/*
 * cubicEquationSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "cubicEquationSolver.hpp"

#include <cmath>

#include "intExpPow.hpp"
#include "exception.hpp"
#include "globalConstants.hpp"
#include "sign.hpp"

std::array<schemi::scalar, 3> schemi::cubicEquationSolver(const scalar Ain,
		const scalar Bin, const scalar Cin, const scalar Din)
{
	typedef long double cScalar;

	constexpr static scalar cubicEquationCriteria { 1E-9 };

	scalar x1 { std::nan("1") }, x2 { std::nan("2") }, x3 { std::nan("3") };

	if ((Ain > cubicEquationCriteria) || (Ain < -cubicEquationCriteria))
	{
		const cScalar a { Bin / Ain }, b { Cin / Ain }, c { Din / Ain };

		const cScalar Q { (pow<cScalar, 2>(a) - 3 * b) / 9 };

		const cScalar R { (a * (2 * pow<cScalar, 2>(a) - 9 * b) + 27 * c) / 54 };

		const cScalar S { pow<cScalar, 3>(Q) - pow<cScalar, 2>(R) };

		if (S >= zeroLevel)
		{
			const cScalar phi { onethirds * std::acos(R / std::sqrt(Q * Q * Q)) };

			x1 = -2. * std::sqrt(Q) * std::cos(phi) - a * onethirds;

			x2 = -2. * std::sqrt(Q) * std::cos(phi - twothirds * Pi_number)
					- a * onethirds;

			x3 = -2. * std::sqrt(Q) * std::cos(phi + twothirds * Pi_number)
					- a * onethirds;
		}
		else if ((S < zeroLevel) && (S >= -zeroLevel))
		{
			x1 = -2. * sign(R) * std::sqrt(Q) - a * onethirds;

			x2 = sign(R) * std::sqrt(Q) - a * onethirds;

			x3 = x2;
		}
		else if (S < -zeroLevel)
			if (Q >= zeroLevel)
			{
				const cScalar phi { onethirds
						* std::acosh(
								std::abs(R)
										/ (std::sqrt(pow<cScalar, 3>(Q))
												+ stabilizator)) };

				x1 = -2 * sign(R) * std::sqrt(Q) * std::cosh(phi)
						- a * onethirds;
			}
			else if ((Q < zeroLevel) && (Q >= -zeroLevel))
				x1 = -std::cbrt(c - pow<cScalar, 3>(a) / 27) - a * onethirds;
			else if (Q < -zeroLevel)
			{
				const cScalar phi { onethirds
						* std::asinh(
								std::abs(R)
										/ (std::sqrt(
												std::abs(pow<cScalar, 3>(Q)))
												+ stabilizator)) };
				x1 = -2 * sign(R) * std::sqrt(std::abs(Q)) * std::sinh(phi)
						- a * onethirds;
			}
			else
				throw exception("Q has unknown value: ",
						errorsEnum::cubicEquationError);
		else
			throw exception("S has unknown value: ",
					errorsEnum::cubicEquationError);
	}
	else
	{
		const cScalar a { Bin }, b { Cin }, c { Din };

		const auto D = pow<cScalar, 2>(b) - 4 * a * c;

		if (D >= 0)
		{
			const auto a2 = 0.5 / (a + stabilizator);

			x1 = a2 * (-b + std::sqrt(D));
			x2 = a2 * (-b - std::sqrt(D));
		}
	}

	return std::array<scalar, 3> { x1, x2, x3 };
}

schemi::scalar schemi::returnSinglePosValue(std::array<scalar, 3> tripleValue)
{
	std::size_t nPositive(0);

	for (auto & v_i : tripleValue)
	{
		if (v_i >= 0)
			nPositive++;
		else
			v_i = -123.;
	}

	if (nPositive == 1)
		return std::max(std::max(tripleValue[0], tripleValue[1]),
				tripleValue[2]);
	else if (nPositive == 2) //Two roots are equal.
	{
		if (((tripleValue[0] - tripleValue[1]) < zeroLevel)
				&& (tripleValue[2] < 0))
			return tripleValue[0];
		else if (((tripleValue[1] - tripleValue[2]) < zeroLevel)
				&& (tripleValue[0] < 0))
			return tripleValue[1];
		else if (((tripleValue[0] - tripleValue[2]) < zeroLevel)
				&& (tripleValue[1] < 0))
			return tripleValue[2];
		else
			throw exception("Can't choose solitary positive value.",
					errorsEnum::positivnessError);
	}
	else if ((nPositive == 3) //All three roots are equal.
			&& (((tripleValue[0] - tripleValue[1]) < zeroLevel)
					&& ((tripleValue[1] - tripleValue[2]) < zeroLevel)))
		return tripleValue[0];
	else
		throw exception("Can't choose solitary positive value.",
				errorsEnum::positivnessError);
}
