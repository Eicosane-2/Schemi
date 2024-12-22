/*
 * biConjugateGradientSolver.hpp
 *
 *  Created on: 2024/12/10
 *      Author: Maxim Boldyrev
 */

#ifndef BICONJUGATEGRADIENTSOLVER_HPP_
#define BICONJUGATEGRADIENTSOLVER_HPP_

#include "abstractMatrixSolver.hpp"

namespace schemi
{
class biConjugateGradientSovler: public abstractMatrixSolver
{
	std::valarray<scalar> algorithm(const std::valarray<scalar> & oldField,
			const SLEMatrix::SLEMatrixStorage & matrix,
			const std::string_view name) const noexcept;
public:
	biConjugateGradientSovler(const std::size_t maxIteration,
			const matrixSolver type_in) noexcept;

	std::valarray<scalar> solve(const std::valarray<scalar> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<vector> solve(const std::valarray<vector> & oldField,
			const SLEMatrix & matrix) const noexcept override;

	std::valarray<tensor> solve(const std::valarray<tensor> & oldField,
			const SLEMatrix & matrix) const noexcept override;
};
}  // namespace schemi

#endif /* BICONJUGATEGRADIENTSOLVER_HPP_ */
