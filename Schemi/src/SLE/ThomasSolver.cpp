/*
 * ThomasSolver.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "ThomasSolver.hpp"

std::valarray<schemi::scalar> schemi::ThomasSolver::algorithm(
		const std::valarray<scalar>&,
		const SLEMatrix::SLEMatrixStorage & matrix,
		const std::string&) const noexcept
{
	std::size_t N { matrix.centralDiagonale.size() };

	std::valarray<scalar> newFieldValue(N);
	std::valarray<scalar> P(N), Q(N);

	P[0] = -matrix.upperTriangle[0][0].first / matrix.centralDiagonale[0];
	Q[0] = matrix.freeTerm[0] / matrix.centralDiagonale[0];
	for (std::size_t i = 1; i < (N - 1); ++i)
	{
		P[i] = -matrix.upperTriangle[i][0].first
				/ (matrix.centralDiagonale[i]
						+ matrix.lowerTriangle[i][0].first * P[i - 1]);
		Q[i] =
				(matrix.freeTerm[i]
						- matrix.lowerTriangle[i][0].first * Q[i - 1])
						/ (matrix.centralDiagonale[i]
								+ matrix.lowerTriangle[i][0].first * P[i - 1]);
	}
	P[N - 1] = 0;
	Q[N - 1] = (matrix.freeTerm[N - 1]
			- matrix.lowerTriangle[N - 1][0].first * Q[N - 2])
			/ (matrix.centralDiagonale[N - 1]
					+ matrix.lowerTriangle[N - 1][0].first * P[N - 2]);

	newFieldValue[N - 1] = Q[N - 1];

	for (std::size_t i = (N - 2);; --i)
	{
		newFieldValue[i] = P[i] * newFieldValue[i + 1] + Q[i];

		if (i == 0)
			break;
	}

	normalize(newFieldValue);

	return newFieldValue;
}

schemi::ThomasSolver::ThomasSolver(const matrixSolver type_in) noexcept :
		abstractMatrixSolver(0, type_in)
{
}

std::valarray<schemi::scalar> schemi::ThomasSolver::solve(
		const std::valarray<scalar> & oldField,
		const SLEMatrix & matrix) const noexcept
{
	return algorithm(oldField, matrix.SLE[0], matrix.name);
}

std::valarray<schemi::vector> schemi::ThomasSolver::solve(
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

std::valarray<schemi::tensor> schemi::ThomasSolver::solve(
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
