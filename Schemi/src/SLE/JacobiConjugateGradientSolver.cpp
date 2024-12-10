/*
 * JacobiConjugateGradientSolver.cpp
 *
 *  Created on: 2023/07/14
 *      Author: Maxim Boldyrev
 */

#include "JacobiConjugateGradientSolver.hpp"

#include <iostream>

#include "SLEMatrixDotProduct.hpp"

std::valarray<schemi::scalar> schemi::JacobiConjugateGradientSolver::algorithm(
		const std::valarray<scalar> & oldField,
		const SLEMatrix::SLEMatrixStorage & matrix,
		const std::string_view name) const noexcept
{
	SLEMatrix::SLEMatrixStorage JacobiPreconditioner = matrix;

	JacobiPreconditioner.centralDiagonale = 1. / matrix.centralDiagonale;

	JacobiPreconditioner.freeTerm = 0.0;

	for (auto & i_row : JacobiPreconditioner.lowerTriangle)
		for (auto & j : i_row)
			j.first = 0;
	for (auto & i_row : JacobiPreconditioner.upperTriangle)
		for (auto & j : i_row)
			j.first = 0;

	std::valarray<scalar> newFieldValues(oldField);

	std::size_t nIterations { 0 };

	std::valarray<scalar> r_n = matrix.freeTerm - (matrix & oldField);
	std::valarray<scalar> d_n = JacobiPreconditioner & r_n;

	while (true)
	{
		nIterations++;

		const scalar diff { std::abs(r_n).max() };

		if ((diff < convergenceTolerance) && (nIterations > 1))
		{
			normalize(newFieldValues);
			return newFieldValues;
		}
		else if (nIterations >= maxIterationNumber)
		{
			std::clog << name << std::endl;
			std::clog
					<< "Conjugate gradient algorithm did not converged. Difference is: "
					<< diff << std::endl;

			normalize(newFieldValues);
			return newFieldValues;
		}
		else
		{
			const scalar alpha = (r_n * (JacobiPreconditioner & r_n)).sum()
					/ ((d_n * (matrix & d_n)).sum() + stabilizator);

			newFieldValues += alpha * d_n;

			const std::valarray<scalar> r_n1 = r_n - alpha * (matrix & d_n);

			const scalar beta =
					(r_n1 * (JacobiPreconditioner & r_n1)).sum()
							/ ((r_n * (JacobiPreconditioner & r_n)).sum()
									+ stabilizator);

			const std::valarray<scalar> d_n1 = (JacobiPreconditioner & r_n1)
					+ beta * d_n;

			r_n = r_n1;
			d_n = d_n1;
		}
	}
}

schemi::JacobiConjugateGradientSolver::JacobiConjugateGradientSolver(
		const std::size_t maxIteration, const matrixSolver type_in) noexcept :
		abstractMatrixSolver(maxIteration, type_in)
{
}

std::valarray<schemi::scalar> schemi::JacobiConjugateGradientSolver::solve(
		const std::valarray<scalar> & oldField,
		const SLEMatrix & matrix) const noexcept
{
	return algorithm(oldField, matrix.SLE[0], matrix.name);
}

std::valarray<schemi::vector> schemi::JacobiConjugateGradientSolver::solve(
		const std::valarray<vector> & oldField,
		const SLEMatrix & matrix) const noexcept
{
	std::valarray<vector> result(oldField.size());

	for (std::size_t j = 0; j < vector::vsize; ++j)
	{
		std::valarray<scalar> v_j_buf(oldField.size());

		for (std::size_t i = 0; i < v_j_buf.size(); ++i)
			v_j_buf[i] = oldField[i]()[j];

		v_j_buf = algorithm(v_j_buf, matrix.SLE[j], matrix.name);

		for (std::size_t i = 0; i < v_j_buf.size(); ++i)
			result[i].r()[j] = v_j_buf[i];
	}

	return result;
}

std::valarray<schemi::tensor> schemi::JacobiConjugateGradientSolver::solve(
		const std::valarray<tensor> & oldField,
		const SLEMatrix & matrix) const noexcept
{
	std::valarray<tensor> result(oldField.size());

	for (std::size_t j = 0; j < tensor::vsize; ++j)
	{
		std::valarray<scalar> v_j_buf(oldField.size());

		for (std::size_t i = 0; i < v_j_buf.size(); ++i)
			v_j_buf[i] = oldField[i]()[j];

		v_j_buf = algorithm(v_j_buf, matrix.SLE[j], matrix.name);

		for (std::size_t i = 0; i < v_j_buf.size(); ++i)
			result[i].r()[j] = v_j_buf[i];
	}

	return result;
}
