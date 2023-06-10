/*
 * noMatrixSolver.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "noMatrixSolver.hpp"

schemi::noMatrixSolver::noMatrixSolver(const matrixSolverEnum type_in) noexcept :
		abstractMatrixSolver(0, type_in)
{
}

std::valarray<schemi::scalar> schemi::noMatrixSolver::solve(
		const std::valarray<scalar> & oldField, const SLEMatrix&) const noexcept
{
	return oldField;
}

std::valarray<schemi::vector> schemi::noMatrixSolver::solve(
		const std::valarray<vector> & oldField, const SLEMatrix&) const noexcept
{
	return oldField;
}

std::valarray<schemi::tensor> schemi::noMatrixSolver::solve(
		const std::valarray<tensor> & oldField, const SLEMatrix&) const noexcept
{
	return oldField;
}
