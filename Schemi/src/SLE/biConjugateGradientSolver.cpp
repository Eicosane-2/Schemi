/*
 * biConjugateGradientSolver.cpp
 *
 *  Created on: 2024/12/10
 *      Author: Maxim Boldyrev
 */

#include "biConjugateGradientSolver.hpp"

#include <iostream>

#include "SLEMatrixDotProduct.hpp"

std::valarray<schemi::scalar> schemi::biConjugateGradientSovler::algorithm(
		const std::valarray<scalar> & oldField,
		const SLEMatrix::SLEMatrixStorage & matrix,
		const std::string_view name) const noexcept
{
	auto matrixT = matrix;
	matrixT.transpose();

	std::valarray<scalar> newFieldValues(oldField);

	std::size_t nIterations { 0 };

	std::valarray<scalar> r_n = matrix.freeTerm - (matrix & oldField);
	std::valarray<scalar> d_n = r_n;

	std::valarray<scalar> ro_n = r_n;
	std::valarray<scalar> do_n = r_n;

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
		//[[unlikely]]
		{
			std::clog << name << std::endl;
			std::clog
					<< "Bi-conjugate gradient algorithm did not converged. Difference is: "
					<< diff << std::endl;

			normalize(newFieldValues);
			return newFieldValues;
		}
		else
		{
			const scalar alpha = (ro_n * r_n).sum()
					/ ((do_n * (matrix & d_n)).sum() + stabilizator);

			newFieldValues += alpha * d_n;

			const std::valarray<scalar> r_n1 = r_n - alpha * (matrix & d_n);
			const std::valarray<scalar> ro_n1 = ro_n - alpha * (matrixT & do_n);

			const scalar beta = (ro_n1 * r_n1).sum()
					/ ((ro_n * r_n).sum() + stabilizator);

			const std::valarray<scalar> d_n1 = r_n1 + beta * d_n;
			const std::valarray<scalar> do_n1 = ro_n1 + beta * do_n;

			r_n = r_n1;
			d_n = d_n1;

			ro_n = ro_n1;
			do_n = do_n1;
		}
	}
}

schemi::biConjugateGradientSovler::biConjugateGradientSovler(
		const std::size_t maxIteration, const matrixSolver type_in) noexcept :
		abstractMatrixSolver(maxIteration, type_in)
{
}

std::valarray<schemi::scalar> schemi::biConjugateGradientSovler::solve(
		const std::valarray<scalar> & oldField,
		const SLEMatrix & matrix) const noexcept
{
	return algorithm(oldField, matrix.SLE[0], matrix.name);
}

std::valarray<schemi::vector> schemi::biConjugateGradientSovler::solve(
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
			result[i].wr()[j] = v_j_buf[i];
	}

	return result;
}

std::valarray<schemi::tensor> schemi::biConjugateGradientSovler::solve(
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
			result[i].wr()[j] = v_j_buf[i];
	}

	return result;
}
