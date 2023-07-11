/*
 * chemicalKineticsChlorumHydrogeniumDissociation.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsChlorumHydrogeniumDissociation.hpp"

#include <iostream>

#include "fieldProducts.hpp"

void schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::normalize(
		std::valarray<scalar> & res) const noexcept
{
	for (auto & i : res)
		if (std::abs(i) < std::numeric_limits<scalar>::epsilon())
			i = 0;
}

schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::cellReactionMatrix() noexcept
{
}

schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::cellReactionMatrix(
		const scalar timeStep, const scalar k_diss_Cl2,
		const scalar k_recomb_Cl2, const scalar k_diss_H2,
		const scalar k_recomb_H2, const scalar C_Cl2_0, const scalar C_Cl_0,
		const scalar C_H2_0, const scalar C_H_0, const scalar M_0,
		const scalar rho_0, const std::valarray<scalar> & molMass) noexcept
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

	matrixDiagonale[0] = A11;
	matrixDiagonale[1] = A22;
	matrixDiagonale[2] = A33;
	matrixDiagonale[3] = A44;

	matrixLeftTriangle[1][0].first = A21;
	matrixLeftTriangle[3][0].first = A43;

	matrixRightTriangle[0][0].first = A12;
	matrixRightTriangle[2][0].first = A34;

	matrixFreeTerm[0] = B1;
	matrixFreeTerm[1] = B2;
	matrixFreeTerm[2] = B3;
	matrixFreeTerm[3] = B4;
}

std::array<schemi::scalar, 4> schemi::chemicalKineticsChlorumHydrogeniumDissociation::cellReactionMatrix::solve(
		const std::array<scalar, 4> & oldField,
		const std::size_t maxIterationNumber) const noexcept
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
				rho, phase.phaseThermodynamics->Mv());
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
		const homogeneousPhase<cubicCell> & phaseIn, const std::size_t itNumber) :
		abstractChemicalKinetics(true, itNumber)
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
