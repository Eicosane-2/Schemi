/*
 * GaussSeidelSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "GaussSeidelSolver.hpp"

#include <iostream>

std::valarray<schemi::scalar> schemi::GaussSeidelSolver::algorithm(
		const std::valarray<scalar> & oldField,
		const SLEMatrix::SLEMatrixStorage & matrix,
		const std::string_view name) const noexcept
{
	std::valarray<scalar> oldIteration(oldField);

	std::valarray<scalar> newIteration(oldIteration);

	std::size_t nIterations { 0 };

	while (true)
	{
		nIterations++;

		for (std::size_t i = 0; i < oldIteration.size(); ++i)
		{
			const scalar aii { 1. / matrix.centralDiagonale[i] };

			const scalar bi { matrix.freeTerm[i] };

			newIteration[i] = bi * aii;

			for (std::size_t j = 0; j < matrix.lowerTriangle[i].size(); ++j)
				newIteration[i] -= matrix.lowerTriangle[i][j].first
						* newIteration[matrix.lowerTriangle[i][j].second] * aii;

			for (std::size_t j = 0; j < matrix.upperTriangle[i].size(); ++j)
				newIteration[i] -= matrix.upperTriangle[i][j].first
						* oldIteration[matrix.upperTriangle[i][j].second] * aii;
		}

		for (std::size_t i = oldIteration.size() - 1;; --i)
		{
			const scalar aii { 1. / matrix.centralDiagonale[i] };

			const scalar bi { matrix.freeTerm[i] };

			newIteration[i] = bi * aii;

			if (matrix.lowerTriangle[i].size() != 0)
				for (std::size_t j = matrix.lowerTriangle[i].size() - 1;; --j)
				{
					newIteration[i] -= matrix.lowerTriangle[i][j].first
							* oldIteration[matrix.lowerTriangle[i][j].second]
							* aii;

					if (j == 0)
						break;
				}

			if (matrix.upperTriangle[i].size() != 0)
				for (std::size_t j = matrix.upperTriangle[i].size() - 1;; --j)
				{
					newIteration[i] -= matrix.upperTriangle[i][j].first
							* newIteration[matrix.upperTriangle[i][j].second]
							* aii;

					if (j == 0)
						break;
				}

			if (i == 0)
				break;
		}

		const scalar diff { relativeIterationDifference(oldIteration,
				newIteration) };

		if (diff < convergenceTolerance)
		{
			normalize(newIteration);
			return newIteration;
		}
		else if (nIterations >= maxIterationNumber)
			[[unlikely]]
			{
				std::clog << name << std::endl;
				std::clog
						<< "Gauss-Seidel algorithm did not converged. Difference is: "
						<< diff << std::endl;

				normalize(newIteration);
				return newIteration;
			}
			else
				oldIteration = newIteration;
		}
	}

	schemi::GaussSeidelSolver::GaussSeidelSolver(const std::size_t maxIteration,
			const matrixSolver type_in) noexcept :
			abstractMatrixSolver(maxIteration, type_in)
	{
	}

	std::valarray<schemi::scalar> schemi::GaussSeidelSolver::solve(
			const std::valarray<scalar> & oldField,
			const SLEMatrix & matrix) const noexcept
	{
		return algorithm(oldField, matrix.SLE[0], matrix.name);
	}

	std::valarray<schemi::vector> schemi::GaussSeidelSolver::solve(
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

	std::valarray<schemi::tensor> schemi::GaussSeidelSolver::solve(
			const std::valarray<tensor> & oldField,
			const SLEMatrix & matrix) const noexcept
	{
		std::valarray<tensor> result(oldField);

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
