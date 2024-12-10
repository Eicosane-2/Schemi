/*
 * abstractChemicalKinetics.hpp
 *
 *  Created on: 2023/05/08
 *      Author: Maxim Boldyrev
 */

#ifndef ABSTRACTCHEMICALKINETICS_HPP_
#define ABSTRACTCHEMICALKINETICS_HPP_

#include <cstddef>
#include <map>
#include <vector>
#include <utility>

#include "chemicalReactionsEnum.hpp"
#include "cubicCell.hpp"
#include "scalar.hpp"
#include "homogeneousPhase.hpp"

namespace schemi
{
namespace chemicalKinetics
{
class abstractChemicalKinetics
{
protected:
	typedef std::vector<std::pair<scalar, std::size_t>> triangleList;

	constexpr static scalar convergenceTolerance { 100
			* convergenceToleranceGlobal };
	constexpr static scalar massFracTolerance { 1E-2 };
	std::size_t maxIterationNumber { 0 };
	const scalar minTimestep { 0 };

	enum class iterativeSolver
	{
		noSolver,
		GaussSeidel,
		ConjugateGradient,
		JacobiConjugateGradient,
		Jacobi,
		GaussElimination
	};

	std::map<std::string, iterativeSolver> solvers;

	struct cellReactingFields
	{
		scalar internalEnergy;

		scalar temperature;

		std::valarray<scalar> concentration;

		std::valarray<scalar> density;

		cellReactingFields(scalar internalEnergy_in, scalar temperature_in,
				std::valarray<scalar> concentration_in,
				std::valarray<scalar> density_in) noexcept :
				internalEnergy(internalEnergy_in), temperature(temperature_in), concentration(
						concentration_in), density(density_in)
		{
		}
	};

	static void normalize(std::valarray<scalar> & res) noexcept;

	template<typename reactionMatrix>
	static auto matrixDotProduct(const reactionMatrix & m,
			const std::valarray<scalar> & v) noexcept->
	std::valarray<scalar>
	{
		std::valarray<scalar> result(v.size());

		for (std::size_t i = 0; i < v.size(); ++i)
		{
			auto i_res = m.Diagonale[i] * v[i];

			const auto & lowTri = m.LeftTriangle[i];

			for (std::size_t j = 0; j < lowTri.size(); ++j)
			{
				const auto index = lowTri[j].second;

				i_res += lowTri[j].first * v[index];
			}

			const auto & upTri = m.RightTriangle[i];

			for (std::size_t j = 0; j < upTri.size(); ++j)
			{
				const auto index = upTri[j].second;

				i_res += upTri[j].first * v[index];
			}

			result[i] = i_res;
		}

		return result;
	}

	template<typename reactionMatrix, std::size_t N>
	static auto solveJ(const reactionMatrix & matrix,
			const std::array<scalar, N> & oldField,
			const std::size_t maxIterationNumber) -> std::array<scalar, N>
	{
		std::valarray<scalar> oldIteration(N);
		for (std::size_t i = 0; i < N; ++i)
			oldIteration[i] = oldField[i];

		std::valarray<scalar> newIteration(oldIteration);

		std::size_t nIterations { 0 };

		while (true)
		{
			nIterations++;

			for (std::size_t i = 0; i < oldIteration.size(); ++i)
			{
				const scalar aii { 1. / matrix.Diagonale[i] };

				const scalar bi { matrix.FreeTerm[i] };

				scalar rowSum { 0 };

				for (std::size_t j = 0; j < matrix.LeftTriangle[i].size(); ++j)
					rowSum -= matrix.LeftTriangle[i][j].first
							* oldIteration[matrix.LeftTriangle[i][j].second];

				for (std::size_t j = 0; j < matrix.RightTriangle[i].size(); ++j)
					rowSum -= matrix.RightTriangle[i][j].first
							* oldIteration[matrix.RightTriangle[i][j].second];

				newIteration[i] = (bi + rowSum) * aii;
			}

			for (std::size_t i = oldIteration.size() - 1;; --i)
			{
				const scalar aii { 1. / matrix.Diagonale[i] };

				const scalar bi { matrix.FreeTerm[i] };

				scalar rowSum { 0 };

				if (matrix.LeftTriangle[i].size() != 0)
					for (std::size_t j = matrix.LeftTriangle[i].size() - 1;;
							--j)
					{
						rowSum -=
								matrix.LeftTriangle[i][j].first
										* oldIteration[matrix.LeftTriangle[i][j].second];

						if (j == 0)
							break;
					}

				if (matrix.RightTriangle[i].size() != 0)
					for (std::size_t j = matrix.RightTriangle[i].size() - 1;;
							--j)
					{
						rowSum -=
								matrix.RightTriangle[i][j].first
										* oldIteration[matrix.RightTriangle[i][j].second];

						if (j == 0)
							break;
					}

				newIteration[i] = (bi + rowSum) * aii;

				if (i == 0)
					break;
			}

			const scalar diff {
					100.
							* std::abs(
									(newIteration - oldIteration)
											/ (std::abs(newIteration)
													+ stabilizator)).max() };

			if (diff < convergenceTolerance)
			{
				normalize(newIteration);

				std::array<scalar, N> ret;
				for (std::size_t i = 0; i < N; ++i)
					ret[i] = newIteration[i];

				return ret;
			}
			else if (nIterations >= maxIterationNumber)
			{
				std::clog
						<< "Jacobi algorithm for chemical reaction did not converged. Difference is: "
						<< diff << std::endl;

				throw exception(
						"Jacobi algorithm for chemical reaction did not converged.",
						errors::systemError);

				normalize(newIteration);

				std::array<scalar, N> ret;
				for (std::size_t i = 0; i < N; ++i)
					ret[i] = newIteration[i];

				return ret;
			}
			else
				oldIteration = newIteration;
		}
	}

	template<typename reactionMatrix, std::size_t N>
	static auto solveGS(const reactionMatrix & matrix,
			const std::array<scalar, N> & oldField,
			const std::size_t maxIterationNumber) -> std::array<scalar, N>
	{
		std::valarray<scalar> oldIteration(N);
		for (std::size_t i = 0; i < N; ++i)
			oldIteration[i] = oldField[i];

		std::valarray<scalar> newIteration(oldIteration);

		std::size_t nIterations { 0 };

		while (true)
		{
			nIterations++;

			for (std::size_t i = 0; i < oldIteration.size(); ++i)
			{
				const scalar aii { 1. / matrix.Diagonale[i] };

				const scalar bi { matrix.FreeTerm[i] };

				newIteration[i] = bi * aii;

				for (std::size_t j = 0; j < matrix.LeftTriangle[i].size(); ++j)
					newIteration[i] -= matrix.LeftTriangle[i][j].first
							* newIteration[matrix.LeftTriangle[i][j].second]
							* aii;

				for (std::size_t j = 0; j < matrix.RightTriangle[i].size(); ++j)
					newIteration[i] -= matrix.RightTriangle[i][j].first
							* oldIteration[matrix.RightTriangle[i][j].second]
							* aii;
			}

			for (std::size_t i = oldIteration.size() - 1;; --i)
			{
				const scalar aii { 1. / matrix.Diagonale[i] };

				const scalar bi { matrix.FreeTerm[i] };

				newIteration[i] = bi * aii;

				if (matrix.LeftTriangle[i].size() != 0)
					for (std::size_t j = matrix.LeftTriangle[i].size() - 1;;
							--j)
					{
						newIteration[i] -= matrix.LeftTriangle[i][j].first
								* oldIteration[matrix.LeftTriangle[i][j].second]
								* aii;

						if (j == 0)
							break;
					}

				if (matrix.RightTriangle[i].size() != 0)
					for (std::size_t j = matrix.RightTriangle[i].size() - 1;;
							--j)
					{
						newIteration[i] -=
								matrix.RightTriangle[i][j].first
										* newIteration[matrix.RightTriangle[i][j].second]
										* aii;

						if (j == 0)
							break;
					}

				if (i == 0)
					break;
			}

			const scalar diff {
					100.
							* std::abs(
									(newIteration - oldIteration)
											/ (std::abs(newIteration)
													+ stabilizator)).max() };

			if (diff < convergenceTolerance)
			{
				normalize(newIteration);

				std::array<scalar, N> ret;
				for (std::size_t i = 0; i < N; ++i)
					ret[i] = newIteration[i];

				return ret;
			}
			else if (nIterations >= maxIterationNumber)
			{
				std::clog
						<< "Gauss-Seidel algorithm for chemical reaction combustion did not converged. Difference is: "
						<< diff << std::endl;

				throw exception(
						"Gauss-Seidel algorithm for chemical reaction combustion did not converged.",
						errors::systemError);

				normalize(newIteration);

				std::array<scalar, N> ret;
				for (std::size_t i = 0; i < N; ++i)
					ret[i] = newIteration[i];

				return ret;
			}
			else
				oldIteration = newIteration;
		}
	}

	template<typename reactionMatrix, std::size_t N>
	static auto solveCG(const reactionMatrix & matrix,
			const std::array<scalar, N> & oldField,
			const std::size_t maxIterationNumber) -> std::array<scalar, N>
	{
		std::valarray<scalar> oldIteration(N);
		for (std::size_t i = 0; i < N; ++i)
			oldIteration[i] = oldField[i];

		std::valarray<scalar> newIteration(oldIteration);

		std::size_t nIterations { 0 };

		std::valarray<scalar> mFrT(N);
		for (std::size_t i = 0; i < N; ++i)
			mFrT[i] = matrix.FreeTerm[i];

		std::valarray<scalar> r_n { mFrT
				- matrixDotProduct(matrix, oldIteration) };
		std::valarray<scalar> d_n = r_n;

		while (true)
		{
			nIterations++;

			const scalar diff { std::abs(r_n).max() };

			if ((diff < convergenceTolerance) && (nIterations > 1))
			{
				normalize(newIteration);

				std::array<scalar, N> ret;
				for (std::size_t i = 0; i < N; ++i)
					ret[i] = newIteration[i];

				return ret;
			}
			else if (nIterations >= maxIterationNumber)
			{
				std::clog
						<< "Conjugate gradient algorithm for chemical reaction combustion did not converged. Difference is: "
						<< diff << std::endl;

				throw exception(
						"Conjugate gradient algorithm for chemical reaction combustion did not converged.",
						errors::systemError);

				normalize(newIteration);

				std::array<scalar, N> ret;
				for (std::size_t i = 0; i < N; ++i)
					ret[i] = newIteration[i];

				return ret;
			}
			else
			{
				const scalar alpha = (r_n * r_n).sum()
						/ ((d_n * matrixDotProduct(matrix, d_n)).sum()
								+ stabilizator);

				newIteration += alpha * d_n;

				const std::valarray<scalar> r_n1 = r_n
						- alpha * matrixDotProduct(matrix, d_n);

				const scalar beta = (r_n1 * r_n1).sum()
						/ ((r_n * r_n).sum() + stabilizator);

				const std::valarray<scalar> d_n1 = r_n1 + beta * d_n;

				r_n = r_n1;
				d_n = d_n1;
			}
		}
	}

	template<typename reactionMatrix, std::size_t N>
	static auto solveJCG(const reactionMatrix & matrix,
			const std::array<scalar, N> & oldField,
			const std::size_t maxIterationNumber) -> std::array<scalar, N>
	{
		reactionMatrix JacobiPreconditioner;

		for (std::size_t i = 0; i < N; ++i)
			JacobiPreconditioner.Diagonale[i] = 1. / matrix.Diagonale[i];

		std::valarray<scalar> oldIteration(N);
		for (std::size_t i = 0; i < N; ++i)
			oldIteration[i] = oldField[i];

		std::valarray<scalar> newIteration(oldIteration);

		std::size_t nIterations { 0 };

		std::valarray<scalar> mFrT(N);
		for (std::size_t i = 0; i < N; ++i)
			mFrT[i] = matrix.FreeTerm[i];

		std::valarray<scalar> r_n { mFrT
				- matrixDotProduct(matrix, oldIteration) };
		std::valarray<scalar> d_n = matrixDotProduct(JacobiPreconditioner, r_n);

		while (true)
		{
			nIterations++;

			const scalar diff { std::abs(r_n).max() };

			if ((diff < convergenceTolerance) && (nIterations > 1))
			{
				normalize(newIteration);

				std::array<scalar, N> ret;
				for (std::size_t i = 0; i < N; ++i)
					ret[i] = newIteration[i];

				return ret;
			}
			else if (nIterations >= maxIterationNumber)
			{
				std::clog
						<< "Jacobi preconditioned conjugate gradient algorithm for chemical reaction combustion did not converged. Difference is: "
						<< diff << std::endl;

				throw exception(
						"Jacobi preconditioned conjugate gradient algorithm for chemical reaction combustion did not converged.",
						errors::systemError);

				normalize(newIteration);

				std::array<scalar, N> ret;
				for (std::size_t i = 0; i < N; ++i)
					ret[i] = newIteration[i];

				return ret;
			}
			else
			{
				oldIteration = newIteration;

				const scalar alpha = (r_n
						* matrixDotProduct(JacobiPreconditioner, r_n)).sum()
						/ ((d_n * matrixDotProduct(matrix, d_n)).sum()
								+ stabilizator);

				newIteration += alpha * d_n;

				const std::valarray<scalar> r_n1 = r_n
						- alpha * matrixDotProduct(matrix, d_n);

				const scalar beta =
						(r_n1 * matrixDotProduct(JacobiPreconditioner, r_n1)).sum()
								/ ((r_n
										* matrixDotProduct(JacobiPreconditioner,
												r_n)).sum() + stabilizator);

				const std::valarray<scalar> d_n1 = matrixDotProduct(
						JacobiPreconditioner, r_n1) + beta * d_n;

				r_n = r_n1;
				d_n = d_n1;
			}
		}
	}

	template<typename reactionMatrix, std::size_t N>
	static auto solveGE(const reactionMatrix & matrix) -> std::array<scalar, N>
	{
		scalar A[N][N];

		auto b = matrix.FreeTerm;

		for (std::size_t i = 0; i < N; ++i)
			A[i][i] = matrix.Diagonale[i];

		for (std::size_t i = 0; i < matrix.LeftTriangle.size(); ++i)
		{
			const auto & lt = matrix.LeftTriangle[i];

			for (std::size_t j = 0; j < lt.size(); ++j)
			{
				const auto & p = lt[j];

				const auto ind = p.second;

				A[i][ind] = p.first;
			}
		}

		for (std::size_t i = 0; i < matrix.RightTriangle.size(); ++i)
		{
			const auto & rt = matrix.RightTriangle[i];

			for (std::size_t j = 0; j < rt.size(); j++)
			{
				const auto & p = rt[j];

				const auto ind = p.second;

				A[i][ind] = p.first;
			}
		}

		for (std::size_t k = 0; k < N - 1; ++k)
			for (std::size_t i = k + 1; i < N; ++i)
			{
				const auto ratio = A[i][k] / A[k][k];

				for (std::size_t j = 0; j < N; ++j)
					A[i][j] = A[i][j] - ratio * A[k][j];

				b[i] = b[i] - ratio * b[k];
			}

		std::valarray<scalar> phi(N);

		phi[N - 1] = b[N - 1] / A[N - 1][N - 1];

		for (std::size_t i = N - 2;; --i)
		{
			scalar term = 0.0;

			for (std::size_t j = i + 1; j < N; ++j)
				term += A[i][j] * phi[j];

			phi[i] = (b[i] - term) / A[i][i];

			if (i == 0)
				break;
		}

		normalize(phi);

		std::array<scalar, N> ret;
		for (std::size_t i = 0; i < N; ++i)
			ret[i] = phi[i];

		return ret;
	}
public:
	const bool chemicalReaction;

	abstractChemicalKinetics(const bool flag, const scalar mt) noexcept;

	virtual ~abstractChemicalKinetics() noexcept =0;

	static std::unique_ptr<abstractChemicalKinetics> createChemicalKinetics(
			const homogeneousPhase<cubicCell> & phaseIn,
			const chemicalReactions chemReactFlag,
			const scalar minimalTimestep) noexcept;

	virtual void solveChemicalKinetics(homogeneousPhase<cubicCell>&) const =0;
};
}
}  // namespace schemi

#endif /* ABSTRACTCHEMICALKINETICS_HPP_ */
