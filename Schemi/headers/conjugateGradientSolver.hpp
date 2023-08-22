/*
 * conjugateGradientSolver.hpp
 *
 *  Created on: 2021/01/23
 *      Author: Maxim Boldyrev
 *
 *      Class for conjugated gradient method for matrix diagonalization.
 */

#ifndef CONJUGATEGRADIENTSOLVER_HPP_
#define CONJUGATEGRADIENTSOLVER_HPP_

#include "abstractMatrixSolver.hpp"

namespace schemi
{
class conjugateGradientSovler: public abstractMatrixSolver
{
	std::valarray<scalar> algorithm(const std::valarray<scalar> & oldField,
			const SLEMatrix::SLEMatrixStorage & matrix,
			const std::string & name) const noexcept;
public:
	conjugateGradientSovler(const std::size_t maxIteration,
			const matrixSolver type_in) noexcept;

	std::valarray<scalar> solve(const std::valarray<scalar> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<vector> solve(const std::valarray<vector> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<tensor> solve(const std::valarray<tensor> & oldField,
			const SLEMatrix & matrix) const noexcept override;
};
}  // namespace schemi

#endif /* CONJUGATEGRADIENTSOLVER_HPP_ */
