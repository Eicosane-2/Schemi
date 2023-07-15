/*
 * chemicalKineticsChlorumHydrogeniumDissociation.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsChlorumHydrogeniumDissociation.hpp"

#include <iostream>

#include "fieldProducts.hpp"

void schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::reactionMatrix::transpose() noexcept
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

void schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::normalize(
		std::valarray<scalar> & res) const noexcept
{
	for (auto & i : res)
		if (std::abs(i) < std::numeric_limits<scalar>::epsilon())
			i = 0;
}

schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::cellReactionMatrix() noexcept :
		solverFlag(iterativeSolver::noSolver), matrix()
{
}

schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::cellReactionMatrix(
		const scalar timeStep, const scalar k_diss_Cl2,
		const scalar k_recomb_Cl2, const scalar k_diss_H2,
		const scalar k_recomb_H2, const scalar C_Cl2_0, const scalar C_Cl_0,
		const scalar C_H2_0, const scalar C_H_0, const scalar M_0,
		const scalar rho_0, const std::valarray<scalar> & molMass,
		const iterativeSolver solverType) :
		solverFlag(solverType), matrix()
{
	const scalar A11 { (1 / timeStep + k_diss_Cl2 * M_0) / molMass[0] };
	const scalar A12 { (-k_recomb_Cl2 * C_Cl_0 * M_0) / molMass[1] };
	const scalar B1 { C_Cl2_0 / (timeStep * rho_0) };

	const scalar A21 { (-2 * k_diss_Cl2 * M_0) / molMass[0] };
	const scalar A22 { (1 / timeStep + 2 * k_recomb_Cl2 * C_Cl_0 * M_0)
			/ molMass[1] };
	const scalar B2 { C_Cl_0 / (timeStep * rho_0) };

	const scalar A33 { (1 / timeStep + k_diss_H2 * M_0) / molMass[2] };
	const scalar A34 { (-k_recomb_H2 * C_H_0 * M_0) / molMass[3] };
	const scalar B3 { C_H2_0 / (timeStep * rho_0) };

	const scalar A43 { (-2 * k_diss_H2 * M_0) / molMass[2] };
	const scalar A44 { (1 / timeStep + 2 * k_recomb_H2 * C_H_0 * M_0)
			/ molMass[3] };
	const scalar B4 { C_H_0 / (timeStep * rho_0) };

	matrix.Diagonale[0] = A11;
	matrix.Diagonale[1] = A22;
	matrix.Diagonale[2] = A33;
	matrix.Diagonale[3] = A44;

	matrix.LeftTriangle[1][0].first = A21;
	matrix.LeftTriangle[3][0].first = A43;

	matrix.RightTriangle[0][0].first = A12;
	matrix.RightTriangle[2][0].first = A34;

	matrix.FreeTerm[0] = B1;
	matrix.FreeTerm[1] = B2;
	matrix.FreeTerm[2] = B3;
	matrix.FreeTerm[3] = B4;

	switch (solverType)
	{
	case iterativeSolver::GaussSeidel:
		break;
	case iterativeSolver::ConjugateGradient:
		break;
	case iterativeSolver::JacobiConjugateGradient:
		break;
	default:
		throw exception(
				"<<solverType>> has been corrupted. Unknown solver type.",
				errorsEnum::initializationError);
		break;
	}
}

std::valarray<schemi::scalar> schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::matrixDotProduct(
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

auto schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::solveGS(
		const std::array<scalar, 4> & oldField,
		const std::size_t maxIterationNumber) const noexcept -> std::array<
scalar, 4>
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
					<< "Gauss-Seidel algorithm for Cl2 and H2 dissociation did not converged. Difference is: "
					<< diff << std::endl;

			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else
			oldIteration = newIteration;
	}
}

auto schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::solveCG(
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
					<< "Conjugate gradient algorithm did not converged for chemical reaction Cl2 and H2 dissociation. Difference is: "
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

auto schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::solveJCG(
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
					<< "Jacobi preconditioned conjugate gradient algorithm did not converged for chemical reaction Cl2 and H2 dissociation. Difference is: "
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
					(rs_n1 * matrixDotProduct(JacobiPreconditioner, rf_n1)).sum()
							/ ((rs_n
									* matrixDotProduct(JacobiPreconditioner,
											rf_n)).sum() + stabilizator);

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

auto schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::solve(
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

std::vector<
		schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix> schemi::chemicalKineticsChlorumHydrogeniumDissociation::velocityCalculation(
		const scalar timestep,
		const homogeneousPhase<cubicCell> & phase) const noexcept
{
	const auto size = phase.pressure.meshRef().cellsSize();

	std::vector<cellReactionMatrix> concentrationVelocityMatrix(size);

	for (std::size_t i = 0; i < size; ++i)
	{
		const scalar T = phase.temperature.ref()[i];

		const scalar k_Cl2_forw = A_Cl2_forw * std::pow(T, n_Cl2_forw)
				* std::exp(-E_Cl2_forw / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_Cl_backw = A_Cl_backw * std::pow(T, n_Cl_backw)
				* std::exp(-E_Cl_backw / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_H2_forw = A_H2_forw * std::pow(T, n_H2_forw)
				* std::exp(-E_H2_forw / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_H_backw = A_H2_backw * std::pow(T, n_H2_backw)
				* std::exp(-E_H2_backw / (phase.phaseThermodynamics->Rv() * T));

		const scalar M = phase.concentration.v[0].ref()[i];
		const scalar Cl2 = phase.concentration.v[1].ref()[i];
		const scalar Cl = phase.concentration.v[2].ref()[i];
		const scalar H2 = phase.concentration.v[3].ref()[i];
		const scalar H = phase.concentration.v[4].ref()[i];

		const scalar rho = phase.density[0].ref()[i];

		concentrationVelocityMatrix[i] = cellReactionMatrix(timestep,
				k_Cl2_forw, k_Cl_backw, k_H2_forw, k_H_backw, Cl2, Cl, H2, H, M,
				rho, phase.phaseThermodynamics->Mv(), itSolv);
	}

	return concentrationVelocityMatrix;
}

void schemi::chemicalKineticsChlorumHydrogeniumDissociation::timeStepIntegration(
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

		auto [newCl2, newCl, newH2, newH] = vC[i].solve(oldMassFraction,
				maxIterationNumber);

		newCl2 = std::max(0., newCl2);
		newCl = std::max(0., newCl);
		newH2 = std::max(0., newH2);
		newH = std::max(0., newH);

		const auto sumFrac = newCl2 + newCl + newH2 + newH;
		newCl2 /= sumFrac;
		newCl /= sumFrac;
		newH2 /= sumFrac;
		newH /= sumFrac;

		const auto newCCl2 = newCl2 * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[0];
		const auto newCCl = newCl * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[1];
		const auto newCH2 = newH2 * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[2];
		const auto newCH = newH * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[3];

		phaseN.concentration.v[1].ref_r()[i] = newCCl2;
		phaseN.concentration.v[2].ref_r()[i] = newCCl;
		phaseN.concentration.v[3].ref_r()[i] = newCH2;
		phaseN.concentration.v[4].ref_r()[i] = newCH;

		phaseN.concentration.v[0].ref_r()[i] = newCCl2 + newCCl + newCH2
				+ newCH;
	}

	phaseN.internalEnergy.ref_r() += deltaU.ref();

	for (std::size_t k = 1; k < phaseN.density.size(); ++k)
		phaseN.density[k].ref_r() = phaseN.concentration.v[k].ref()
				* phaseN.phaseThermodynamics->Mv()[k - 1];
}

schemi::chemicalKineticsChlorumHydrogeniumDissociation::chemicalKineticsChlorumHydrogeniumDissociation(
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

	chem >> skipBuffer >> A_Cl2_forw;
	chem >> skipBuffer >> n_Cl2_forw;
	chem >> skipBuffer >> E_Cl2_forw;

	chem >> skipBuffer >> A_Cl_backw;
	chem >> skipBuffer >> n_Cl_backw;
	chem >> skipBuffer >> E_Cl_backw;

	chem >> skipBuffer >> A_H2_forw;
	chem >> skipBuffer >> n_H2_forw;
	chem >> skipBuffer >> E_H2_forw;

	chem >> skipBuffer >> A_H2_backw;
	chem >> skipBuffer >> n_H2_backw;
	chem >> skipBuffer >> E_H2_backw;

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
}

void schemi::chemicalKineticsChlorumHydrogeniumDissociation::solveChemicalKinetics(
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
