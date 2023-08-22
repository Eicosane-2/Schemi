/*
 * chemicalKineticsNO2Disproportionation.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsNO2Disproportionation.hpp"

#include <iostream>
#include <numeric>

void schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::reactionMatrix::transpose() noexcept
{
	std::array<triangleList, 4> LeftTriangleNew, RightTriangleNew;

	for (std::size_t i = 0; i < Diagonale.size(); ++i)
	{
		for (std::size_t j = 0; j < LeftTriangle[i].size(); ++j)
		{
			const std::size_t jAbsOld = LeftTriangle[i][j].second;
			const auto Avalue = LeftTriangle[i][j].first;

			const std::size_t iNew = jAbsOld;
			const std::size_t jNew = i;

			RightTriangleNew[iNew].push_back( { Avalue, jNew });
		}

		for (std::size_t j = 0; j < RightTriangle[i].size(); ++j)
		{
			const std::size_t jAbsOld = RightTriangle[i][j].second;
			const auto Avalue = RightTriangle[i][j].first;

			const std::size_t iNew = jAbsOld;
			const std::size_t jNew = i;

			LeftTriangleNew[iNew].push_back( { Avalue, jNew });
		}
	}

	LeftTriangle = LeftTriangleNew;
	RightTriangle = RightTriangleNew;
}

void schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::normalize(
		std::valarray<scalar> & res) const noexcept
{
	for (auto & i : res)
		if (std::abs(i) < std::numeric_limits<scalar>::epsilon())
			i = 0;
}

schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::cellReactionMatrix() noexcept :
		solverFlag(iterativeSolver::noSolver), matrix()
{
}

schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::cellReactionMatrix(
		const scalar timeStep, const scalar k_f, const scalar k_b,
		const scalar C_NO2_0, const scalar C_H2O_0, const scalar C_HNO2_0,
		const scalar C_HNO3_0, const scalar rho_0,
		const std::array<scalar, 4> & molMass, const iterativeSolver solverType) :
		solverFlag(solverType), matrix()
{
	const scalar A11 { (1 / timeStep + k_f * C_NO2_0 * C_H2O_0) / molMass[0] };

	const scalar A12 { (k_f * C_NO2_0 * C_NO2_0) / molMass[1] };
	const scalar A13 { -(k_b * C_HNO3_0) / molMass[2] };
	const scalar A14 { -(k_b * C_HNO2_0) / molMass[3] };

	const scalar B1 { C_NO2_0 / (timeStep * rho_0) };

	const scalar A21 { (0.5 * k_f * C_NO2_0 * C_H2O_0) / molMass[0] };

	const scalar A22 { (1 / timeStep + 0.5 * k_f * C_NO2_0 * C_NO2_0)
			/ molMass[1] };

	const scalar A23 { (-0.5 * k_b * C_HNO3_0) / molMass[2] };
	const scalar A24 { (-0.5 * k_b * C_HNO2_0) / molMass[3] };

	const scalar B2 { C_H2O_0 / (timeStep * rho_0) };

	const scalar A31 { (-0.5 * k_f * C_NO2_0 * C_H2O_0) / molMass[0] };
	const scalar A32 { (-0.5 * k_f * C_NO2_0 * C_NO2_0) / molMass[1] };

	const scalar A33 { (1 / timeStep + 0.5 * k_b * C_HNO3_0) / molMass[2] };

	const scalar A34 { (0.5 * k_b * C_HNO2_0) / molMass[3] };

	const scalar B3 { C_HNO2_0 / (timeStep * rho_0) };

	const scalar A41 { (-0.5 * k_f * C_NO2_0 * C_H2O_0) / molMass[0] };
	const scalar A42 { (-0.5 * k_f * C_NO2_0 * C_NO2_0) / molMass[1] };
	const scalar A43 { (0.5 * k_b * C_HNO3_0) / molMass[2] };
	const scalar A44 { (1 / timeStep + 0.5 * k_b * C_HNO2_0) / molMass[3] };

	const scalar B4 { C_HNO3_0 / (timeStep * rho_0) };

	matrix.Diagonale[0] = A11;
	matrix.Diagonale[1] = A22;
	matrix.Diagonale[2] = A33;
	matrix.Diagonale[3] = A44;

	matrix.LeftTriangle[1][0].first = A21;

	matrix.LeftTriangle[2][0].first = A31;
	matrix.LeftTriangle[2][1].first = A32;

	matrix.LeftTriangle[3][0].first = A41;
	matrix.LeftTriangle[3][1].first = A42;
	matrix.LeftTriangle[3][2].first = A43;

	matrix.RightTriangle[0][0].first = A12;
	matrix.RightTriangle[0][1].first = A13;
	matrix.RightTriangle[0][2].first = A14;

	matrix.RightTriangle[1][0].first = A23;
	matrix.RightTriangle[1][1].first = A24;

	matrix.RightTriangle[2][0].first = A34;

	matrix.FreeTerm[0] = B1;
	matrix.FreeTerm[1] = B2;
	matrix.FreeTerm[2] = B3;
	matrix.FreeTerm[3] = B4;
}

std::valarray<schemi::scalar> schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::matrixDotProduct(
		const reactionMatrix & m,
		const std::valarray<scalar> & v) const noexcept
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

auto schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::solveJ(
		const std::array<scalar, 4> & oldField,
		const std::size_t maxIterationNumber) const -> std::array<
		schemi::scalar, 4>
{
	std::valarray<scalar> oldIteration { oldField[0], oldField[1], oldField[2],
			oldField[3] };

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
				for (std::size_t j = matrix.LeftTriangle[i].size() - 1;; --j)
				{
					rowSum -= matrix.LeftTriangle[i][j].first
							* oldIteration[matrix.LeftTriangle[i][j].second];

					if (j == 0)
						break;
				}

			if (matrix.RightTriangle[i].size() != 0)
				for (std::size_t j = matrix.RightTriangle[i].size() - 1;; --j)
				{
					rowSum -= matrix.RightTriangle[i][j].first
							* oldIteration[matrix.RightTriangle[i][j].second];

					if (j == 0)
						break;
				}

			newIteration[i] = (bi + rowSum) * aii;

			if (i == 0)
				break;
		}

		const scalar diff { 100.
				* std::abs(
						(newIteration - oldIteration)
								/ (std::abs(newIteration) + stabilizator)).max() };

		if (diff < convergenceTolerance)
		{
			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else if (nIterations >= maxIterationNumber)
		{
			std::clog
					<< "Jacobi algorithm for NO2 disproportionation did not converged. Difference is: "
					<< diff << std::endl;

			throw exception(
					"Jacobi algorithm for NO2 disproportionation did not converged.",
					errors::systemError);

			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else
			oldIteration = newIteration;
	}
}

auto schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::solveGS(
		const std::array<scalar, 4> & oldField,
		const std::size_t maxIterationNumber) const -> std::array<
		schemi::scalar, 4>
{
	std::valarray<scalar> oldIteration { oldField[0], oldField[1], oldField[2],
			oldField[3] };

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
						* newIteration[matrix.LeftTriangle[i][j].second] * aii;

			for (std::size_t j = 0; j < matrix.RightTriangle[i].size(); ++j)
				newIteration[i] -= matrix.RightTriangle[i][j].first
						* oldIteration[matrix.RightTriangle[i][j].second] * aii;
		}

		for (std::size_t i = oldIteration.size() - 1;; --i)
		{
			const scalar aii { 1. / matrix.Diagonale[i] };

			const scalar bi { matrix.FreeTerm[i] };

			newIteration[i] = bi * aii;

			if (matrix.LeftTriangle[i].size() != 0)
				for (std::size_t j = matrix.LeftTriangle[i].size() - 1;; --j)
				{
					newIteration[i] -= matrix.LeftTriangle[i][j].first
							* oldIteration[matrix.LeftTriangle[i][j].second]
							* aii;

					if (j == 0)
						break;
				}

			if (matrix.RightTriangle[i].size() != 0)
				for (std::size_t j = matrix.RightTriangle[i].size() - 1;; --j)
				{
					newIteration[i] -= matrix.RightTriangle[i][j].first
							* newIteration[matrix.RightTriangle[i][j].second]
							* aii;

					if (j == 0)
						break;
				}

			if (i == 0)
				break;
		}

		const scalar diff { 100.
				* std::abs(
						(newIteration - oldIteration)
								/ (std::abs(newIteration) + stabilizator)).max() };

		if (diff < convergenceTolerance)
		{
			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else if (nIterations >= maxIterationNumber)
		{
			std::clog
					<< "Gauss-Seidel algorithm for NO2 disproportionation did not converged. Difference is: "
					<< diff << std::endl;

			throw exception(
					"Gauss-Seidel algorithm for NO2 disproportionation did not converged.",
					errors::systemError);

			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else
			oldIteration = newIteration;
	}
}

auto schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::solveCG(
		const std::array<scalar, 4> & oldField,
		const std::size_t maxIterationNumber) const -> std::array<
		scalar, 4>
{
	std::valarray<scalar> oldIteration { oldField[0], oldField[1], oldField[2],
			oldField[3] };

	std::valarray<scalar> newIteration(oldIteration);

	auto matrixT = matrix;
	matrixT.transpose();

	std::size_t nIterations { 0 };

	std::valarray<scalar> rf_n { std::valarray<scalar> { matrix.FreeTerm[0],
			matrix.FreeTerm[1], matrix.FreeTerm[2], matrix.FreeTerm[3] }
			- matrixDotProduct(matrix, oldIteration) };
	std::valarray<scalar> df_n = rf_n;

	auto rs_n = rf_n;
	auto ds_n = df_n;

	while (true)
	{
		nIterations++;

		//const scalar diff { 100.
		//		* std::abs(
		//				(newIteration - oldIteration)
		//						/ (std::abs(newIteration) + stabilizator)).max() };

		const scalar diff { rf_n.max() };

		if ((diff < convergenceTolerance) && (nIterations > 1))
		{
			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else if (nIterations >= maxIterationNumber)
		{
			std::clog
					<< "Conjugate gradient algorithm did not converged for chemical reaction NO2 disproportionation. Difference is: "
					<< diff << std::endl;

			throw exception(
					"Conjugate gradient algorithm for NO2 disproportionation did not converged.",
					errors::systemError);

			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else
		{
			oldIteration = newIteration;

			const scalar alpha = (rs_n * rf_n).sum()
					/ ((ds_n * matrixDotProduct(matrix, df_n)).sum()
							+ stabilizator);

			newIteration += alpha * df_n;

			const std::valarray<scalar> rf_n1 = rf_n
					- alpha * matrixDotProduct(matrix, df_n);
			const std::valarray<scalar> rs_n1 = rs_n
					- alpha * matrixDotProduct(matrixT, ds_n);

			const scalar beta = (rs_n1 * rf_n1).sum()
					/ ((rs_n * rf_n).sum() + stabilizator);

			const std::valarray<scalar> df_n1 = rf_n1 + beta * df_n;
			const std::valarray<scalar> ds_n1 = rs_n1 + beta * ds_n;

			rf_n = rf_n1;
			df_n = df_n1;

			rs_n = rs_n1;
			ds_n = ds_n1;
		}
	}
}

auto schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::solveJCG(
		const std::array<scalar, 4> & oldField,
		const std::size_t maxIterationNumber) const -> std::array<
		scalar, 4>
{
	reactionMatrix JacobiPreconditioner;

	JacobiPreconditioner.Diagonale = { 1. / matrix.Diagonale[0], 1.
			/ matrix.Diagonale[1], 1. / matrix.Diagonale[2], 1.
			/ matrix.Diagonale[3] };

	std::valarray<scalar> oldIteration { oldField[0], oldField[1], oldField[2],
			oldField[3] };

	std::valarray<scalar> newIteration(oldIteration);

	auto matrixT = matrix;
	matrixT.transpose();

	std::size_t nIterations { 0 };

	std::valarray<scalar> rf_n { std::valarray<scalar> { matrix.FreeTerm[0],
			matrix.FreeTerm[1], matrix.FreeTerm[2], matrix.FreeTerm[3] }
			- matrixDotProduct(matrix, oldIteration) };
	std::valarray<scalar> df_n = matrixDotProduct(JacobiPreconditioner, rf_n);

	auto rs_n = rf_n;
	auto ds_n = matrixDotProduct(JacobiPreconditioner, rs_n);

	while (true)
	{
		nIterations++;

		//const scalar diff { 100.
		//		* std::abs(
		//				(newIteration - oldIteration)
		//						/ (std::abs(newIteration) + stabilizator)).max() };

		const scalar diff { rf_n.max() };

		if ((diff < convergenceTolerance) && (nIterations > 1))
		{
			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else if (nIterations >= maxIterationNumber)
		{
			std::clog
					<< "Jacobi preconditioned conjugate gradient algorithm did not converged for chemical reaction NO2 disproportionation. Difference is: "
					<< diff << std::endl;

			throw exception(
					"Jacobi preconditioned conjugate gradient algorithm for NO2 disproportionation did not converged.",
					errors::systemError);

			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else
		{
			oldIteration = newIteration;

			const scalar alpha = (rs_n
					* matrixDotProduct(JacobiPreconditioner, rf_n)).sum()
					/ ((ds_n * matrixDotProduct(matrix, df_n)).sum()
							+ stabilizator);

			newIteration += alpha * df_n;

			const std::valarray<scalar> rf_n1 = rf_n
					- alpha * matrixDotProduct(matrix, df_n);
			const std::valarray<scalar> rs_n1 = rs_n
					- alpha * matrixDotProduct(matrixT, ds_n);

			const scalar beta =
					(rs_n1 * matrixDotProduct(JacobiPreconditioner, rf_n1)).sum()
							/ ((rs_n
									* matrixDotProduct(JacobiPreconditioner,
											rs_n)).sum() + stabilizator);

			const std::valarray<scalar> df_n1 = matrixDotProduct(
					JacobiPreconditioner, rf_n1) + beta * df_n;
			const std::valarray<scalar> ds_n1 = matrixDotProduct(
					JacobiPreconditioner, rs_n1) + beta * ds_n;

			rf_n = rf_n1;
			df_n = df_n1;

			rs_n = rs_n1;
			ds_n = ds_n1;
		}
	}
}

auto schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::solveGE() const ->
std::array<scalar, 4>
{
	constexpr std::size_t N { 4 };

	scalar A[N][N] { { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0,
			0 } };

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

	std::valarray<scalar> phi(4);

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

	return
	{	phi[0], phi[1], phi[2], phi[3]};
}

auto schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::solve(
		const std::array<scalar, 4> & oldField,
		const std::size_t maxIterationNumber) const -> std::array<
		scalar, 4>
{
	switch (solverFlag)
	{
	case iterativeSolver::GaussSeidel:
		return solveGS(oldField, maxIterationNumber);
		break;
	case iterativeSolver::ConjugateGradient:
		return solveCG(oldField, maxIterationNumber);
		break;
	case iterativeSolver::JacobiConjugateGradient:
		return solveJCG(oldField, maxIterationNumber);
		break;
	case iterativeSolver::Jacobi:
		return solveJ(oldField, maxIterationNumber);
		break;
	case iterativeSolver::GaussElimination:
		return solveGE();
		break;
	default:
		throw exception("Unknown chemical iterative solver type.",
				errors::initializationError);
		break;
	}
}

schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix schemi::chemicalKineticsNO2Disproportionation::velocityCalculation(
		const scalar timestep, const scalar T,
		const std::array<scalar, 5> & concentrations,
		const std::array<scalar, 4> & molarMasses, const scalar rho,
		const scalar R) const noexcept
{
	const scalar k_forward = A_forward * std::pow(T, n_forward)
			* std::exp(-E_forward / (R * T));

	const scalar k_backward = A_backward * std::pow(T, n_backward)
			* std::exp(-E_backward / (R * T));

	const scalar & NO2 = concentrations[1];
	const scalar & H2O = concentrations[2];
	const scalar & HNO2 = concentrations[3];
	const scalar & HNO3 = concentrations[4];

	return cellReactionMatrix(timestep, k_forward, k_backward, NO2, H2O, HNO2,
			HNO3, rho, molarMasses, itSolv);
}

void schemi::chemicalKineticsNO2Disproportionation::timeStepIntegration(
		homogeneousPhase<cubicCell> & phaseN) const
{
	auto & mesh_ = phaseN.pressure.meshRef();

	const scalar timeStep = mesh_.timestep();

	std::size_t subItNum = 1;

	auto subTimeStep = timeStep / subItNum;

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		std::valarray<scalar> concs(phaseN.concentration.v.size()), dens(
				phaseN.density.size() - 1);

		for (std::size_t k = 0; k < concs.size(); ++k)
			concs[k] = phaseN.concentration.v[k].ref()[i];

		for (std::size_t k = 0; k < dens.size(); ++k)
			dens[k] = phaseN.density[k + 1].ref()[i];

		const cellReactingFields oldValues(phaseN.internalEnergy.ref()[i],
				phaseN.temperature.ref()[i], concs, dens);

		cellReactingFields newValues(oldValues);

		while (true)
		{
			try
			{
				newValues = oldValues;

				for (std::size_t st = 0; st < subItNum; ++st)
				{
					scalar deltaU { 0 };

					const std::array<scalar, 4> oldMassFraction {
							newValues.density[0] / phaseN.density[0].ref()[i],
							newValues.density[1] / phaseN.density[0].ref()[i],
							newValues.density[2] / phaseN.density[0].ref()[i],
							newValues.density[3] / phaseN.density[0].ref()[i] };

					const scalar sumFracOld = std::accumulate(
							oldMassFraction.begin(), oldMassFraction.end(), 0.);

					if (sumFracOld == 0.0)
						continue;

					const auto cellReactionVel = velocityCalculation(
							subTimeStep, phaseN.temperature.ref()[i],
							{ newValues.concentration[0],
									newValues.concentration[1],
									newValues.concentration[2],
									newValues.concentration[3],
									newValues.concentration[4] },
							{ phaseN.phaseThermodynamics->Mv()[0],
									phaseN.phaseThermodynamics->Mv()[1],
									phaseN.phaseThermodynamics->Mv()[2],
									phaseN.phaseThermodynamics->Mv()[3] },
							phaseN.density[0].ref()[i],
							phaseN.phaseThermodynamics->Rv());

					auto [newNO2, newH2O, newHNO2, newHNO3] =
							cellReactionVel.solve(oldMassFraction,
									maxIterationNumber);

					newNO2 = std::max(static_cast<scalar>(0.), newNO2);
					newH2O = std::max(static_cast<scalar>(0.), newH2O);
					newHNO2 = std::max(static_cast<scalar>(0.), newHNO2);
					newHNO3 = std::max(static_cast<scalar>(0.), newHNO3);

					const auto sumFracNew = newNO2 + newH2O + newHNO2 + newHNO3;

					if (std::abs(sumFracNew - sumFracOld) > massFracTolerance)
						throw exception(
								std::string(
										"Error in mass fraction is too big. Delta is ")
										+ std::to_string(
												std::abs(
														sumFracNew
																- sumFracOld))
										+ '.', errors::systemError);

					newNO2 /= sumFracNew;
					newH2O /= sumFracNew;
					newHNO2 /= sumFracNew;
					newHNO3 /= sumFracNew;

					newNO2 *= sumFracOld;
					newH2O *= sumFracOld;
					newHNO2 *= sumFracOld;
					newHNO3 *= sumFracOld;

					const auto newCNO2 = newNO2 * phaseN.density[0].ref()[i]
							/ phaseN.phaseThermodynamics->Mv()[0];
					const auto newCH2O = newH2O * phaseN.density[0].ref()[i]
							/ phaseN.phaseThermodynamics->Mv()[1];
					const auto newCHNO2 = newHNO2 * phaseN.density[0].ref()[i]
							/ phaseN.phaseThermodynamics->Mv()[2];
					const auto newCHNO3 = newHNO3 * phaseN.density[0].ref()[i]
							/ phaseN.phaseThermodynamics->Mv()[3];

					{
						const scalar deltaC_HNO3 = newCHNO3
								- phaseN.concentration.v[4].ref_r()[i];

						const auto & thermo = *phaseN.phaseThermodynamics;

						const auto deltaCv = thermo.Cvv()[3] + thermo.Cvv()[2]
								- 2 * thermo.Cvv()[0] - thermo.Cvv()[1];

						deltaU = -(ΔU_298
								+ deltaCv * (newValues.temperature - 298.15))
								* deltaC_HNO3;
					}

					newValues.concentration[1] = newCNO2;
					newValues.concentration[2] = newCH2O;
					newValues.concentration[3] = newCHNO2;
					newValues.concentration[4] = newCHNO3;

					newValues.concentration[0] = 0;
					for (std::size_t k = 1; k < newValues.concentration.size();
							++k)
						newValues.concentration[0] +=
								newValues.concentration[k];

					newValues.internalEnergy += deltaU;

					newValues.temperature = phaseN.phaseThermodynamics->TFromUv(
							newValues.concentration, newValues.internalEnergy);

					if (newValues.temperature < 0)
						throw exception(
								"Negative temperature after chemical reaction.",
								errors::negativeTemperatureError);

					for (std::size_t k = 0; k < newValues.density.size(); ++k)
						newValues.density[k] = newValues.concentration[k + 1]
								* phaseN.phaseThermodynamics->Mv()[k];
				}
				break;
			} catch (const exception & ex)
			{
				subItNum *= 10;
				subTimeStep = timeStep / subItNum;

				std::cout << ex.what() << std::endl;
				std::cout
						<< "Restarting chemical reaction step with lower timestep. Timestep decreased in "
						<< subItNum << " times." << std::endl;

				if (subTimeStep < minTimestep)
					throw exception("Timestep become too small.",
							errors::systemError);
			}
		}

		phaseN.internalEnergy.ref_r()[i] = newValues.internalEnergy;
		phaseN.temperature.ref_r()[i] = newValues.temperature;
		for (std::size_t k = 0; k < newValues.concentration.size(); ++k)
			phaseN.concentration.v[k].ref_r()[i] = newValues.concentration[k];
		for (std::size_t k = 0; k < newValues.density.size(); ++k)
			phaseN.density[k + 1].ref_r()[i] = newValues.density[k];
	}
}

schemi::chemicalKineticsNO2Disproportionation::chemicalKineticsNO2Disproportionation(
		const homogeneousPhase<cubicCell> & phaseIn, const scalar mt) :
		abstractChemicalKinetics(true, mt), itSolv(iterativeSolver::noSolver)
{
	if (phaseIn.concentration.v.size() < 5)
		throw exception("Wrong number of substances.",
				errors::initializationError);

	std::string skipBuffer;

	std::ifstream chem { "./set/chemicalKinetics.txt" };

	if (chem.is_open())
		std::cout << "./set/chemicalKinetics.txt is opened." << std::endl;
	else
		throw exception("./set/chemicalKinetics.txt not found.",
				errors::initializationError);

	chem >> skipBuffer >> skipBuffer;

	chem >> skipBuffer >> A_forward;
	chem >> skipBuffer >> n_forward;
	chem >> skipBuffer >> E_forward;

	chem >> skipBuffer >> A_backward;
	chem >> skipBuffer >> n_backward;
	chem >> skipBuffer >> E_backward;

	chem >> skipBuffer >> ΔН_298;

	std::string solverName;

	chem >> skipBuffer >> solverName;

	if (solverName == "Gauss-Seidel")
		itSolv = iterativeSolver::GaussSeidel;
	else if (solverName == "Conjugate_gradient")
		itSolv = iterativeSolver::ConjugateGradient;
	else if (solverName == "Jacobi_gonjugate_gradient")
		itSolv = iterativeSolver::JacobiConjugateGradient;
	else if (solverName == "Jacobi")
		itSolv = iterativeSolver::Jacobi;
	else if (solverName == "Gauss_elimination")
		itSolv = iterativeSolver::GaussElimination;
	else
		throw exception("Unknown type of chemical iterative solver.",
				errors::initializationError);

	chem >> skipBuffer >> maxIterationNumber;

	ΔU_298 = ΔН_298 - Δn * phaseIn.phaseThermodynamics->Rv() * 298.15;
}

void schemi::chemicalKineticsNO2Disproportionation::solveChemicalKinetics(
		homogeneousPhase<cubicCell> & phaseIn) const
{
	auto phaseN1 = phaseIn;

	timeStepIntegration(phaseN1);

	timeStepIntegration(phaseN1);

	phaseN1.pressure.ref_r() = phaseN1.phaseThermodynamics->pFromUv(
			phaseN1.concentration.p, phaseN1.internalEnergy.ref());

	phaseN1.HelmholtzEnergy.ref_r() = phaseN1.phaseThermodynamics->Fv(
			phaseN1.concentration.p, phaseN1.temperature.ref());

	phaseN1.entropy.ref_r() = phaseN1.phaseThermodynamics->Sv(
			phaseN1.concentration.p, phaseN1.temperature.ref());

	{
		const auto v2 = ampProduct(phaseN1.velocity, phaseN1.velocity);

		phaseN1.totalEnergy.ref_r() = phaseN1.internalEnergy.ref()
				+ phaseN1.density[0].ref() * v2.ref() * 0.5
				+ phaseN1.rhokTurb.ref();
	}

	phaseIn.average(phaseN1, *phaseIn.phaseThermodynamics, 0.5);
}
