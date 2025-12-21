/*
 * abstractMatrixSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractMatrixSolver.hpp"

#include <algorithm>
#include <cmath>

#include "biConjugateGradientSolver.hpp"
#include "conjugateGradientSolver.hpp"
#include "GaussSeidelSolver.hpp"
#include "globalConstants.hpp"
#include "JacobiConjugateGradientSolver.hpp"
#include "JacobiSolver.hpp"
#include "noMatrixSolver.hpp"
#include "ThomasSolver.hpp"

schemi::scalar schemi::abstractMatrixSolver::relativeIterationDifference(
		const std::valarray<scalar> & oldField,
		const std::valarray<scalar> & newField) const noexcept
{
	return 100.
			* std::abs(
					(newField - oldField) / (std::abs(newField) + stabilizator)).max();
}

void schemi::abstractMatrixSolver::normalize(
		std::valarray<scalar> & res) const noexcept
{
	std::replace_if(std::begin(res), std::end(res), [](const auto & i) 
	{
		return std::abs(i) < std::numeric_limits<scalar>::epsilon();
	}, 0);
}

schemi::abstractMatrixSolver::abstractMatrixSolver(const std::size_t maxIt_in,
		const matrixSolver type_in) noexcept :
		maxIterationNumber(maxIt_in), solverType(type_in)
{
}

schemi::abstractMatrixSolver::~abstractMatrixSolver() noexcept
{
}

std::pair<std::unique_ptr<schemi::abstractMatrixSolver>,
		std::unique_ptr<schemi::abstractMatrixSolver>> schemi::abstractMatrixSolver::createMatrixSolver(
		const std::string_view name, const std::string_view dim,
		const std::size_t iter)
{
	matrixSolver matrixSolverFlag;
	if (name == "GaussSeidel")
		matrixSolverFlag = matrixSolver::GaussSeidel;
	else if ((name == "Thomas") && (dim == "1D"))
		matrixSolverFlag = matrixSolver::Thomas;
	else if (name == "noSolver")
		matrixSolverFlag = matrixSolver::noSovler;
	else if (name == "explicit")
		matrixSolverFlag = matrixSolver::explicitSolver;
	else if (name == "conjugateGradient")
		matrixSolverFlag = matrixSolver::conjugateGradient;
	else if (name == "JacobiConjugateGradient")
		matrixSolverFlag = matrixSolver::JacobiConjugateGradient;
	else if (name == "Jacobi")
		matrixSolverFlag = matrixSolver::Jacobi;
	else if (name == "biConjugateGradient")
		matrixSolverFlag = matrixSolver::biConjugateGradient;
	else
		[[unlikely]]
		throw exception("Unknown or inappropriate matrix solver flag.",
				errors::initialisationError);

	switch (matrixSolverFlag)
	{
	case matrixSolver::conjugateGradient:
		return std::make_pair(
				std::make_unique<conjugateGradientSovler>(iter,
						matrixSolverFlag),
				std::make_unique<conjugateGradientSovler>(iter,
						matrixSolverFlag));
		break;
	case matrixSolver::JacobiConjugateGradient:
		return std::make_pair(
				std::make_unique<JacobiConjugateGradientSolver>(iter,
						matrixSolverFlag),
				std::make_unique<JacobiConjugateGradientSolver>(iter,
						matrixSolverFlag));
		break;
	case matrixSolver::Thomas:
		return std::make_pair(std::make_unique<ThomasSolver>(matrixSolverFlag),
				std::make_unique<GaussSeidelSolver>(iter,
						matrixSolver::GaussSeidel));
		break;
	case matrixSolver::explicitSolver:
	case matrixSolver::noSovler:
		return std::make_pair(
				std::make_unique<noMatrixSolver>(matrixSolverFlag),
				std::make_unique<noMatrixSolver>(matrixSolverFlag));
		break;
	case matrixSolver::Jacobi:
		return std::make_pair(
				std::make_unique<JacobiSolver>(iter, matrixSolver::Jacobi),
				std::make_unique<JacobiSolver>(iter, matrixSolver::Jacobi));
		break;
	case matrixSolver::biConjugateGradient:
		return std::make_pair(
				std::make_unique<biConjugateGradientSovler>(iter,
						matrixSolver::biConjugateGradient),
				std::make_unique<biConjugateGradientSovler>(iter,
						matrixSolver::biConjugateGradient));
		break;
	default:
		return std::make_pair(
				std::make_unique<GaussSeidelSolver>(iter,
						matrixSolver::GaussSeidel),
				std::make_unique<GaussSeidelSolver>(iter,
						matrixSolver::GaussSeidel));
		break;
	}
}
