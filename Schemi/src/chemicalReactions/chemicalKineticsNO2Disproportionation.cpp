/*
 * chemicalKineticsNO2Disproportionation.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsNO2Disproportionation.hpp"

#include <iostream>

void schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::normalize(
		std::valarray<scalar> & res) const noexcept
{
	for (auto & i : res)
		if (std::abs(i) < std::numeric_limits<scalar>::epsilon())
			i = 0;
}

schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::cellReactionMatrix() noexcept
{
}

schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::cellReactionMatrix(
		const scalar timeStep, const scalar k_f, const scalar k_b,
		const scalar C_NO2_0, const scalar C_H2O_0, const scalar C_HNO2_0,
		const scalar C_HNO3_0, const scalar rho_0,
		const std::valarray<scalar> & molMass) noexcept
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

	matrixDiagonale[0] = A11;
	matrixDiagonale[1] = A22;
	matrixDiagonale[2] = A33;
	matrixDiagonale[3] = A44;

	matrixLeftTriangle[1][0].first = A21;

	matrixLeftTriangle[2][0].first = A31;
	matrixLeftTriangle[2][1].first = A32;

	matrixLeftTriangle[3][0].first = A41;
	matrixLeftTriangle[3][1].first = A42;
	matrixLeftTriangle[3][2].first = A43;

	matrixRightTriangle[0][0].first = A12;
	matrixRightTriangle[0][1].first = A13;
	matrixRightTriangle[0][2].first = A14;

	matrixRightTriangle[1][0].first = A23;
	matrixRightTriangle[1][1].first = A24;

	matrixRightTriangle[2][0].first = A34;

	matrixFreeTerm[0] = B1;
	matrixFreeTerm[1] = B2;
	matrixFreeTerm[2] = B3;
	matrixFreeTerm[3] = B4;
}

std::array<schemi::scalar, 4> schemi::chemicalKineticsNO2Disproportionation::cellReactionMatrix::solve(
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
					<< "Gauss-Seidel algorithm for H2 + Cl2 combustion did not converged. Difference is: "
					<< diff << std::endl;

			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3]};
		}
		else
			oldIteration = newIteration;
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
				phase.phaseThermodynamics->Mv());
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

	chem >> skipBuffer >> A_forward;
	chem >> skipBuffer >> n_forward;
	chem >> skipBuffer >> E_forward;

	chem >> skipBuffer >> A_backward;
	chem >> skipBuffer >> n_backward;
	chem >> skipBuffer >> E_backward;

	chem >> skipBuffer >> ΔН_298;

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
