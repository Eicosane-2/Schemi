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
	explicit noMatrixSolver(const matrixSolver type_in) noexcept;

	std::valarray<scalar> solve(const std::valarray<scalar> & oldField,
			const SLEMatrix&) const noexcept override;

	std::valarray<vector> solve(const std::valarray<vector> & oldField,
			const SLEMatrix&) const noexcept override;

	std::valarray<tensor> solve(const std::valarray<tensor> & oldField,
			const SLEMatrix&) const noexcept override;
};
}  // namespace schemi

#endif /* NOMATRIXSOLVER_HPP_ */
