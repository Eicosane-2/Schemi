/*
 * JacobiSolver.hpp
 *
 *  Created on: 2023/07/15
 *      Author: Maxim Boldyrev
 */

#ifndef JACOBISOLVER_HPP_
#define JACOBISOLVER_HPP_

#include "abstractMatrixSolver.hpp"

namespace schemi
{
class JacobiSolver: public abstractMatrixSolver
{
	std::valarray<scalar> algorithm(const std::valarray<scalar> & oldField,
			const SLEMatrix::SLEMatrixStorage & matrix,
			const std::string & name) const noexcept;
public:
	JacobiSolver(const std::size_t maxIteration,
			const matrixSolverEnum type_in) noexcept;

	std::valarray<scalar> solve(const std::valarray<scalar> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<vector> solve(const std::valarray<vector> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<tensor> solve(const std::valarray<tensor> & oldField,
			const SLEMatrix & matrix) const noexcept override;
};
}  // namespace schemi

#endif /* JACOBISOLVER_HPP_ */
