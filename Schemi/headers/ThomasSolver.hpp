/*
 * ThomasSolver.hpp
 *
 *  Created on: 2020/02/22
 *      Author: Maxim Boldyrev
 *
 *      Class for Thomas solver for matrix of system of linear equations (only for tridiagonal).
 */

#ifndef THOMASSOLVER_HPP_
#define THOMASSOLVER_HPP_

#include "abstractMatrixSolver.hpp"

namespace schemi
{
class ThomasSolver: public abstractMatrixSolver
{
	std::valarray<scalar> algorithm(const std::valarray<scalar>&,
			const SLEMatrix::SLEMatrixStorage & matrix,
			const std::string&) const noexcept;
public:
	explicit ThomasSolver(const matrixSolverEnum type_in) noexcept;

	std::valarray<scalar> solve(const std::valarray<scalar> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<vector> solve(const std::valarray<vector> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<tensor> solve(const std::valarray<tensor> & oldField,
			const SLEMatrix & matrix) const noexcept override;
};
}  // namespace schemi

#endif /* THOMASSOLVER_HPP_ */
