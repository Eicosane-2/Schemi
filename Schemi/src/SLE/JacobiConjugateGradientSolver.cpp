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
		const std::string & name) const noexcept
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

	std::valarray<scalar> oldFieldValues(oldField);

	std::valarray<scalar> newFieldValues(oldFieldValues);

	auto matrixT = matrix;
	matrixT.transpose();

	std::size_t nIterations { 0 };

	std::valarray<scalar> rf_n = matrix.freeTerm - (matrix & oldFieldValues);
	std::valarray<scalar> df_n = JacobiPreconditioner & rf_n;

	auto rs_n = rf_n;
	auto ds_n = JacobiPreconditioner & rs_n;

	while (true)
	{
		nIterations++;

		//const scalar diff { relativeIterationDifference(oldFieldValues,
		//		newFieldValues) };

		const scalar diff { rf_n.max() };

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
			oldFieldValues = newFieldValues;

			const scalar alpha = (rs_n * (JacobiPreconditioner & rf_n)).sum()
					/ ((ds_n * (matrix & df_n)).sum() + stabilizator);

			newFieldValues += alpha * df_n;

			const std::valarray<scalar> rf_n1 = rf_n - alpha * (matrix & df_n);
			const std::valarray<scalar> rs_n1 = rs_n - alpha * (matrixT & ds_n);

			const scalar beta = (rs_n1 * (JacobiPreconditioner & rf_n1)).sum()
					/ ((rs_n * (JacobiPreconditioner & rf_n)).sum()
							+ stabilizator);

			const std::valarray<scalar> df_n1 = (JacobiPreconditioner & rf_n1)
					+ beta * df_n;
			const std::valarray<scalar> ds_n1 = (JacobiPreconditioner & rs_n1)
					+ beta * ds_n;

			rf_n = rf_n1;
			df_n = df_n1;

			rs_n = rs_n1;
			ds_n = ds_n1;
		}
	}
}

schemi::JacobiConjugateGradientSolver::JacobiConjugateGradientSolver(
		const std::size_t maxIteration, const matrixSolverEnum type_in) noexcept :
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
			v_j_buf[i] = oldField[i].v()[j];

		v_j_buf = algorithm(v_j_buf, matrix.SLE[j], matrix.name);

		for (std::size_t i = 0; i < v_j_buf.size(); ++i)
			result[i].v_r()[j] = v_j_buf[i];
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
			v_j_buf[i] = oldField[i].v()[j];

		v_j_buf = algorithm(v_j_buf, matrix.SLE[j], matrix.name);

		for (std::size_t i = 0; i < v_j_buf.size(); ++i)
			result[i].v_r()[j] = v_j_buf[i];
	}

	return result;
}
