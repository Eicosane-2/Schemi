/*
 * GaussSeidelSolver.hpp
 *
 *  Created on: 2020/01/19
 *      Author: Maxim Boldyrev
 *
 *      Class for Gauss-Seidel solver for matrix of system of linear equations.
 */

#ifndef GAUSSSEIDELSOLVER_HPP_
#define GAUSSSEIDELSOLVER_HPP_

#include "abstractMatrixSolver.hpp"

namespace schemi
{
class GaussSeidelSolver: public abstractMatrixSolver
{
	std::valarray<scalar> algorithm(const std::valarray<scalar> & oldField,
			const SLEMatrix::SLEMatrixStorage & matrix,
			const std::string & name) const noexcept;
public:
	GaussSeidelSolver(const std::size_t maxIteration,
			const matrixSolverEnum type_in) noexcept;

	std::valarray<scalar> solve(const std::valarray<scalar> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<vector> solve(const std::valarray<vector> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<tensor> solve(const std::valarray<tensor> & oldField,
			const SLEMatrix & matrix) const noexcept override;
};
}  // namespace schemi

#endif /* GAUSSSEIDELSOLVER_HPP_ */
