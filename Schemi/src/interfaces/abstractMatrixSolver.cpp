/*
 * abstractMatrixSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractMatrixSolver.hpp"

#include <cmath>

#include "globalConstants.hpp"

schemi::scalar schemi::abstractMatrixSolver::relativeIterationDifference(
		const std::valarray<scalar> & oldField,
		const std::valarray<scalar> & newField) const noexcept
{
	return 100.
			* std::abs(
					(newField - oldField) / (std::abs(newField) + stabilizator)).max();
}

void schemi::abstractMatrixSolver::normalize(
		std::valarray<scalar> & res) const noexcept
{
	for (auto & r_i : res)
		if (std::abs(r_i) < std::numeric_limits<scalar>::epsilon())
			r_i = 0;
}

schemi::abstractMatrixSolver::abstractMatrixSolver(const std::size_t maxIt_in,
		const matrixSolver type_in) noexcept :
		maxIterationNumber(maxIt_in), solverType(type_in)
{
}

schemi::abstractMatrixSolver::~abstractMatrixSolver() noexcept
{
}
