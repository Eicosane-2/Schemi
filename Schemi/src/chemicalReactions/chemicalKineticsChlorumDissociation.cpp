/*
 * chemicalKineticsChlorumDissociation.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsChlorumDissociation.hpp"

#include <iostream>

void schemi::chemicalKineticsChlorumDissociation::cellReactionMatrix::normalize(
		std::valarray<scalar> & res) const noexcept
{
	for (auto & i : res)
		if (std::abs(i) < std::numeric_limits<scalar>::epsilon())
			i = 0;
}

schemi::chemicalKineticsChlorumDissociation::cellReactionMatrix::cellReactionMatrix() noexcept
{
}

schemi::chemicalKineticsChlorumDissociation::cellReactionMatrix::cellReactionMatrix(
		const scalar timeStep, const scalar k_diss, const scalar k_recomb,
		const scalar C_Cl2_0, const scalar C_Cl_0, const scalar M_0,
		const scalar rho_0, const std::valarray<scalar> & molMass) noexcept
{
	const scalar A11 { (1 / timeStep + k_diss * M_0) / molMass[0] };
	const scalar A12 { (-k_recomb * C_Cl_0 * M_0) / molMass[1] };
	const scalar B1 { C_Cl2_0 / (rho_0 * timeStep) };

	const scalar A21 { (-2 * k_diss * M_0) / molMass[0] };
	const scalar A22 { (1 / timeStep + 2 * k_recomb * C_Cl_0 * M_0) / molMass[1] };
	const scalar B2 { C_Cl_0 / (timeStep * rho_0) };

	matrixDiagonale[0] = A11;
	matrixDiagonale[1] = A22;

	matrixRightTriangle[0][0].first = A12;
	matrixLeftTriangle[1][0].first = A21;

	matrixFreeTerm[0] = B1;
	matrixFreeTerm[1] = B2;
}

std::array<schemi::scalar, 2> schemi::chemicalKineticsChlorumDissociation::cellReactionMatrix::solve(
		const std::array<scalar, 2> & oldField,
		const std::size_t maxIterationNumber) const noexcept
{
	std::valarray<scalar> oldIteration { oldField[0], oldField[1] };

	std::valarray<scalar> newIteration(oldIteration);

	std::size_t nIterations { 0 };

	while (true)
	{
		nIterations++;

		for (std::size_t i = 0; i < oldIteration.size(); ++i)
		{
			const scalar aii { 1. / matrixDiagonale[i] };

			const scalar bi { matrixFreeTerm[i] };

			newIteration[i] = bi * aii;

			for (std::size_t j = 0; j < matrixLeftTriangle[i].size(); ++j)
				newIteration[i] -= matrixLeftTriangle[i][j].first
						* newIteration[matrixLeftTriangle[i][j].second] * aii;

			for (std::size_t j = 0; j < matrixRightTriangle[i].size(); ++j)
				newIteration[i] -= matrixRightTriangle[i][j].first
						* oldIteration[matrixRightTriangle[i][j].second] * aii;
		}

		for (std::size_t i = oldIteration.size() - 1;; --i)
		{
			const scalar aii { 1. / matrixDiagonale[i] };

			const scalar bi { matrixFreeTerm[i] };

			newIteration[i] = bi * aii;

			if (matrixLeftTriangle[i].size() != 0)
				for (std::size_t j = matrixLeftTriangle[i].size() - 1;; --j)
				{
					newIteration[i] -= matrixLeftTriangle[i][j].first
							* oldIteration[matrixLeftTriangle[i][j].second]
							* aii;

					if (j == 0)
						break;
				}

			if (matrixRightTriangle[i].size() != 0)
				for (std::size_t j = matrixRightTriangle[i].size() - 1;; --j)
				{
					newIteration[i] -= matrixRightTriangle[i][j].first
							* newIteration[matrixRightTriangle[i][j].second]
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
			{	newIteration[0], newIteration[1]};
		}
		else if (nIterations >= maxIterationNumber)
		{
			std::clog
					<< "Gauss-Seidel algorithm for Cl2 dissociation did not converged. Difference is: "
					<< diff << std::endl;

			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1]};
		}
		else
			oldIteration = newIteration;
	}
}

std::vector<schemi::chemicalKineticsChlorumDissociation::cellReactionMatrix> schemi::chemicalKineticsChlorumDissociation::velocityCalculation(
		const scalar timestep,
		const homogeneousPhase<cubicCell> & phase) const noexcept
{
	const auto size = phase.pressure.meshRef().cellsSize();

	std::vector<cellReactionMatrix> concentrationVelocityMatrix(size);

	for (std::size_t i = 0; i < size; ++i)
	{
		const scalar T = phase.temperature.ref()[i];

		const scalar k_forw = A_forw * std::pow(T, n_forw)
				* std::exp(-E_forw / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_backw = A_backw * std::pow(T, n_backw)
				* std::exp(-E_backw / (phase.phaseThermodynamics->Rv() * T));

		const scalar Cl2 = phase.concentration.v[1].ref()[i];
		const scalar Cl = phase.concentration.v[2].ref()[i];
		const scalar M = phase.concentration.v[0].ref()[i];

		const scalar rho = phase.density[0].ref()[i];

		concentrationVelocityMatrix[i] = cellReactionMatrix(timestep, k_forw,
				k_backw, Cl2, Cl, M, rho, phase.phaseThermodynamics->Mv());
	}

	return concentrationVelocityMatrix;
}

void schemi::chemicalKineticsChlorumDissociation::timeStepIntegration(
		homogeneousPhase<cubicCell> & phaseN) const noexcept
{
	auto & mesh_ = phaseN.pressure.meshRef();

	const scalar timeStep = mesh_.timestep();

	volumeField<scalar> deltaU(mesh_, 0);

	const auto vC = velocityCalculation(timeStep, phaseN);

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const std::array<scalar, 2> oldMassFraction {
				phaseN.concentration.v[1].ref()[i]
						* phaseN.phaseThermodynamics->Mv()[0]
						/ phaseN.density[0].ref()[i],
				phaseN.concentration.v[2].ref()[i]
						* phaseN.phaseThermodynamics->Mv()[1]
						/ phaseN.density[0].ref()[i] };

		auto [newCl2, newCl] = vC[i].solve(oldMassFraction, maxIterationNumber);

		newCl2 = std::max(0., newCl2);
		newCl = std::max(0., newCl);

		const auto sumFrac = newCl2 + newCl;
		newCl2 /= sumFrac;
		newCl /= sumFrac;

		const auto newCCl2 = newCl2 * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[0];
		const auto newCCl = newCl * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[1];

		phaseN.concentration.v[1].ref_r()[i] = newCCl2;
		phaseN.concentration.v[2].ref_r()[i] = newCCl;

		phaseN.concentration.v[0].ref_r()[i] = newCCl2 + newCCl;
	}

	phaseN.internalEnergy.ref_r() += deltaU.ref();

	for (std::size_t k = 1; k < phaseN.density.size(); ++k)
		phaseN.density[k].ref_r() = phaseN.concentration.v[k].ref()
				* phaseN.phaseThermodynamics->Mv()[k - 1];
}

schemi::chemicalKineticsChlorumDissociation::chemicalKineticsChlorumDissociation(
		const homogeneousPhase<cubicCell> & phaseIn, const std::size_t itNumber) :
		abstractChemicalKinetics(true, itNumber)
{
	if (phaseIn.concentration.v.size() < 3)
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

	chem >> skipBuffer >> A_forw;
	chem >> skipBuffer >> n_forw;
	chem >> skipBuffer >> E_forw;

	chem >> skipBuffer >> A_backw;
	chem >> skipBuffer >> n_backw;
	chem >> skipBuffer >> E_backw;
}

void schemi::chemicalKineticsChlorumDissociation::solveChemicalKinetics(
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
