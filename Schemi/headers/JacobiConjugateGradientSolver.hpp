/*
 * JacobiConjugateGradientSolver.hpp
 *
 *  Created on: 2023/07/14
 *      Author: Maxim Boldyrev
 *
 *      Class for conjugated gradient with Jacobi preconditioner method for matrix diagonalization.
 */

#ifndef JACOBICONJUGATEGRADIENTSOLVER_HPP_
#define JACOBICONJUGATEGRADIENTSOLVER_HPP_

#include "abstractMatrixSolver.hpp"

namespace schemi
{
class JacobiConjugateGradientSolver: public abstractMatrixSolver
{
	std::valarray<scalar> algorithm(const std::valarray<scalar> & oldField,
			const SLEMatrix::SLEMatrixStorage & matrix,
			const std::string & name) const noexcept;
public:
	JacobiConjugateGradientSolver(const std::size_t maxIteration,
			const matrixSolverEnum type_in) noexcept;

	std::valarray<scalar> solve(const std::valarray<scalar> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<vector> solve(const std::valarray<vector> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<tensor> solve(const std::valarray<tensor> & oldField,
			const SLEMatrix & matrix) const noexcept override;
};
}  // namespace schemi

#endif /* JACOBICONJUGATEGRADIENTSOLVER_HPP_ */
