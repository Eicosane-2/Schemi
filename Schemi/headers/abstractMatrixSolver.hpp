/*
 * abstractMatrixSolver.hpp
 *
 *  Created on: 2020/02/21
 *      Author: Maxim Boldyrev
 *
 *      Interface class for solver of system of linear equations' matrix.
 */

#ifndef ABSTRACTMATRIXSOLVER_HPP_
#define ABSTRACTMATRIXSOLVER_HPP_

#include <cstddef>
#include <valarray>

#include "matrixSolverEnum.hpp"

#include "scalar.hpp"
#include "SLEMatrix.hpp"
#include "vector.hpp"
#include "tensor.hpp"

namespace schemi
{
class abstractMatrixSolver
{
protected:
	constexpr static scalar convergenceTolerance { convergenceToleranceGlobal };
	const std::size_t maxIterationNumber;

	scalar relativeIterationDifference(const std::valarray<scalar> & oldField,
			const std::valarray<scalar> & newField) const noexcept;

	void normalize(std::valarray<scalar> & res) const noexcept;
public:
	virtual ~abstractMatrixSolver() noexcept =0;

	static std::pair<std::unique_ptr<abstractMatrixSolver>,
			std::unique_ptr<abstractMatrixSolver>> createMatrixSolver(
			const std::string_view name, const std::string_view dim,
			const std::size_t iter);

	abstractMatrixSolver(const std::size_t maxIt_in,
			const matrixSolver type_in) noexcept;

	virtual std::valarray<scalar> solve(
			const std::valarray<scalar>& /*oldField*/,
			const SLEMatrix& /*matrix*/) const noexcept =0;

	virtual std::valarray<vector> solve(
			const std::valarray<vector>& /*oldField*/,
			const SLEMatrix& /*matrix*/) const noexcept =0;

	virtual std::valarray<tensor> solve(
			const std::valarray<tensor>& /*oldField*/,
			const SLEMatrix& /*matrix*/) const noexcept =0;

	const matrixSolver solverType;
};
}  // namespace schemi

#endif /* ABSTRACTMATRIXSOLVER_HPP_ */
