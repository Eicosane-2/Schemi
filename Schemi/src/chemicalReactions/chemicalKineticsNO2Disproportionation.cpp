/*
 * chemicalKineticsNO2Disproportionation.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsNO2Disproportionation.hpp"

#include <iostream>

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
		const std::valarray<scalar> & molMass, const iterativeSolver solverType) :
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

auto schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::solveGS(
		const std::array<scalar, 4> & oldField,
		const std::size_t maxIterationNumber) const noexcept -> std::array<
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
		const std::size_t maxIterationNumber) const noexcept -> std::array<
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

		const scalar diff { 100.
				* std::abs(
						(newIteration - oldIteration)
								/ (std::abs(newIteration) + stabilizator)).max() };

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
		const std::size_t maxIterationNumber) const noexcept -> std::array<
scalar, 4>
{
	reactionMatrix JacobiPreconditioner;

	JacobiPreconditioner.Diagonale = { matrix.Diagonale[0], matrix.Diagonale[1],
			matrix.Diagonale[2], matrix.Diagonale[3] };

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

		const scalar diff { 100.
				* std::abs(
						(newIteration - oldIteration)
								/ (std::abs(newIteration) + stabilizator)).max() };

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
					(rf_n1 * matrixDotProduct(JacobiPreconditioner, rf_n1)).sum()
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
	default:
		throw exception("Unknown chemical iterative solver type.",
				errorsEnum::initializationError);
		break;
	}
}

std::vector<schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix> schemi::chemicalKineticsNO2Disproportionation::velocityCalculation(
		const scalar timestep,
		const homogeneousPhase<cubicCell> & phase) const noexcept
{
	const auto size = phase.pressure.meshRef().cellsSize();

	std::vector<cellReactionMatrix> concentrationVelocityMatrix(size);

	for (std::size_t i = 0; i < size; ++i)
	{
		const scalar T = phase.temperature.ref()[i];

		const scalar k_forward = A_forward * std::pow(T, n_forward)
				* std::exp(-E_forward / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_backward = A_backward * std::pow(T, n_backward)
				* std::exp(-E_backward / (phase.phaseThermodynamics->Rv() * T));

		const scalar NO2 = phase.concentration.v[1].ref()[i];
		const scalar H2O = phase.concentration.v[2].ref()[i];
		const scalar HNO2 = phase.concentration.v[3].ref()[i];
		const scalar HNO3 = phase.concentration.v[4].ref()[i];

		const scalar rho = phase.density[0].ref()[i];

		concentrationVelocityMatrix[i] = cellReactionMatrix(timestep, k_forward,
				k_backward, NO2, H2O, HNO2, HNO3, rho,
				phase.phaseThermodynamics->Mv(), itSolv);
	}

	return concentrationVelocityMatrix;
}

void schemi::chemicalKineticsNO2Disproportionation::timeStepIntegration(
		homogeneousPhase<cubicCell> & phaseN) const noexcept
{
	auto & mesh_ = phaseN.pressure.meshRef();

	const scalar timeStep = mesh_.timestep();

	volumeField<scalar> deltaU(mesh_, 0);

	std::valarray<scalar> sourceTimeStep(veryBig, mesh_.cellsSize());

	const auto vC = velocityCalculation(timeStep, phaseN);

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const std::array<scalar, 4> oldMassFraction {
				phaseN.concentration.v[1].ref()[i]
						* phaseN.phaseThermodynamics->Mv()[0]
						/ phaseN.density[0].ref()[i],
				phaseN.concentration.v[2].ref()[i]
						* phaseN.phaseThermodynamics->Mv()[1]
						/ phaseN.density[0].ref()[i],
				phaseN.concentration.v[3].ref()[i]
						* phaseN.phaseThermodynamics->Mv()[2]
						/ phaseN.density[0].ref()[i],
				phaseN.concentration.v[4].ref()[i]
						* phaseN.phaseThermodynamics->Mv()[3]
						/ phaseN.density[0].ref()[i] };

		auto [newNO2, newH2O, newHNO2, newHNO3] = vC[i].solve(oldMassFraction,
				maxIterationNumber);

		newNO2 = std::max(0., newNO2);
		newH2O = std::max(0., newH2O);
		newHNO2 = std::max(0., newHNO2);
		newHNO3 = std::max(0., newHNO3);

		const auto sumFrac = newNO2 + newH2O + newHNO2 + newHNO3;
		newNO2 /= sumFrac;
		newH2O /= sumFrac;
		newHNO2 /= sumFrac;
		newHNO3 /= sumFrac;

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

			deltaU.ref_r()[i] = -(ΔU_298
					+ deltaCv * (phaseN.temperature.ref()[i] - 298.15))
					* deltaC_HNO3;
		}

		phaseN.concentration.v[1].ref_r()[i] = newCNO2;
		phaseN.concentration.v[2].ref_r()[i] = newCH2O;
		phaseN.concentration.v[3].ref_r()[i] = newCHNO2;
		phaseN.concentration.v[4].ref_r()[i] = newCHNO3;

		phaseN.concentration.v[0].ref_r()[i] = newCNO2 + newCH2O + newCHNO2
				+ newCHNO3;
	}

	phaseN.internalEnergy.ref_r() += deltaU.ref();

	for (std::size_t k = 1; k < phaseN.density.size(); ++k)
		phaseN.density[k].ref_r() = phaseN.concentration.v[k].ref()
				* phaseN.phaseThermodynamics->Mv()[k - 1];
}

schemi::chemicalKineticsNO2Disproportionation::chemicalKineticsNO2Disproportionation(
		const homogeneousPhase<cubicCell> & phaseIn) :
		abstractChemicalKinetics(true), itSolv(iterativeSolver::noSolver)
{
	if (phaseIn.concentration.v.size() < 5)
		throw exception("Wrong number of substances.",
				errorsEnum::initializationError);

	std::string skipBuffer;

	std::ifstream chem { "./set/chemicalKinetics.txt" };

	if (chem.is_open())
		std::cout << "./set/chemicalKinetics.txt is opened." << std::endl;
	else
		throw exception("./set/chemicalKinetics.txt not found.",
				errorsEnum::initializationError);

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
	else if (solverName == "Conjugate_Gradient")
		itSolv = iterativeSolver::ConjugateGradient;
	else if (solverName == "Jacobi_Conjugate_Gradient")
		itSolv = iterativeSolver::JacobiConjugateGradient;
	else
		throw exception("Unknown type of chemical iterative solver.",
				errorsEnum::initializationError);

	chem >> skipBuffer >> maxIterationNumber;

	ΔU_298 = ΔН_298 - Δn * phaseIn.phaseThermodynamics->Rv() * 298.15;
}

void schemi::chemicalKineticsNO2Disproportionation::solveChemicalKinetics(
		homogeneousPhase<cubicCell> & phaseIn) const noexcept
{
	auto phaseN1 = phaseIn;

	timeStepIntegration(phaseN1);

	phaseN1.temperature.ref_r() = phaseN1.phaseThermodynamics->TFromUv(
			phaseN1.concentration.p, phaseN1.internalEnergy.ref());

	timeStepIntegration(phaseN1);

	phaseN1.temperature.ref_r() = phaseN1.phaseThermodynamics->TFromUv(
			phaseN1.concentration.p, phaseN1.internalEnergy.ref());

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
