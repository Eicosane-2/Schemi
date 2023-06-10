/*
 * noMatrixSolver.hpp
 *
 *  Created on: 2020/02/29
 *      Author: Maxim Boldyrev
 *
 *      Class for not diagonalization of matrix.
 */

#ifndef NOMATRIXSOLVER_HPP_
#define NOMATRIXSOLVER_HPP_

#include "abstractMatrixSolver.hpp"

namespace schemi
{
class noMatrixSolver: public abstractMatrixSolver
{
public:
	explicit noMatrixSolver(const matrixSolverEnum type_in) noexcept;

	std::valarray<scalar> solve(const std::valarray<scalar> & oldField,
			const SLEMatrix&) const noexcept;

	std::valarray<vector> solve(const std::valarray<vector> & oldField,
			const SLEMatrix&) const noexcept;

	std::valarray<tensor> solve(const std::valarray<tensor> & oldField,
			const SLEMatrix&) const noexcept;
};
}  // namespace schemi

#endif /* NOMATRIXSOLVER_HPP_ */
