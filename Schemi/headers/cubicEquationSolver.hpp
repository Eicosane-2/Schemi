/*
 * cubicEquationSolver.hpp
 *
 *  Created on: 2020/02/23
 *      Author: Maxim Boldyrev
 *
 *      Functions for cubic equation solver: trigonometric solver and single positive value finder.
 */

#ifndef CUBICEQUATIONSOLVER_HPP_
#define CUBICEQUATIONSOLVER_HPP_

#include <array>

#include "scalar.hpp"

namespace schemi
{
/*Viete's trigonometric formula*/
std::array<scalar, 3> cubicEquationSolver(const scalar Ain, const scalar Bin,
		const scalar Cin, const scalar Din);

scalar returnSinglePosValue(std::array<scalar, 3> tripleValue);
}  // namespace schemi

#endif /* CUBICEQUATIONSOLVER_HPP_ */
