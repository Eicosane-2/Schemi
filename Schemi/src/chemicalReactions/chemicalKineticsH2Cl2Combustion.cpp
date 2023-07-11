/*
 * chemicalKineticsH2Cl2Combustion.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsH2Cl2Combustion.hpp"

#include <iostream>

void schemi::chemicalKineticsH2Cl2Combustion::cellReactionMatrix::normalize(
		std::valarray<scalar> & res) const noexcept
{
	for (auto & i : res)
		if (std::abs(i) < std::numeric_limits<scalar>::epsilon())
			i = 0;
}

schemi::chemicalKineticsH2Cl2Combustion::cellReactionMatrix::cellReactionMatrix() noexcept
{
}

schemi::chemicalKineticsH2Cl2Combustion::cellReactionMatrix::cellReactionMatrix(
		const scalar timeStep, const scalar k_diss_Cl2,
		const scalar k_recomb_Cl, const scalar k_diss_H2,
		const scalar k_recomb_H, const scalar k_prop_Cl_H2,
		const scalar k_prop_H_Cl2, const scalar k_recomb_H_Cl,
		const scalar k_diss_HCl, const scalar C_Cl2_0, const scalar C_Cl_0,
		const scalar C_H2_0, const scalar C_H_0, const scalar C_HCl_0,
		const scalar M_0, const scalar rho_0,
		const std::valarray<scalar> & molMass) noexcept
{
	const scalar A11 { (1 / timeStep + k_diss_Cl2 * M_0
			+ 0.5 * k_prop_H_Cl2 * C_H_0) / molMass[0] };

	const scalar A12 { (-k_recomb_Cl * C_Cl_0 * M_0) / molMass[1] };
	const scalar A14 { (0.5 * k_prop_H_Cl2 * C_Cl2_0) / molMass[3] };

	const scalar B1 { C_Cl2_0 / (timeStep * rho_0) };

	const scalar A21 { (-2 * k_diss_Cl2 * M_0 - 0.5 * k_prop_H_Cl2 * C_H_0)
			/ molMass[0] };

	const scalar A22 { (1 / timeStep + 2 * k_recomb_Cl * C_Cl_0 * M_0
			+ 0.5 * k_prop_Cl_H2 * C_H2_0 + 0.5 * k_recomb_H_Cl * C_H_0 * M_0)
			/ molMass[1] };

	const scalar A23 { (0.5 * k_prop_Cl_H2 * C_Cl_0) / molMass[2] };
	const scalar A24 { (-0.5 * k_prop_H_Cl2 * C_Cl2_0
			+ 0.5 * k_recomb_H_Cl * C_Cl_0 * M_0) / molMass[3] };
	const scalar A25 { (-k_diss_HCl * M_0) / molMass[4] };

	const scalar B2 { C_Cl_0 / (timeStep * rho_0) };

	const scalar A32 { (0.5 * k_prop_Cl_H2 * C_H2_0) / molMass[1] };

	const scalar A33 { (1 / timeStep + k_diss_H2 * M_0
			+ 0.5 * k_prop_Cl_H2 * C_Cl_0) / molMass[2] };

	const scalar A34 { (-k_recomb_H * C_H_0 * M_0) / molMass[3] };

	const scalar B3 { C_H2_0 / (timeStep * rho_0) };

	const scalar A41 { (0.5 * k_prop_H_Cl2 * C_H_0) / molMass[0] };
	const scalar A42 { (-0.5 * k_prop_Cl_H2 * C_H2_0
			+ 0.5 * k_recomb_H_Cl * C_H_0 * M_0) / molMass[1] };
	const scalar A43 { (-2 * k_diss_H2 * M_0 - 0.5 * k_prop_Cl_H2 * C_Cl_0)
			/ molMass[2] };

	const scalar A44 { (1 / timeStep + 2 * k_recomb_H * C_H_0 * M_0
			+ 0.5 * k_prop_H_Cl2 * C_Cl2_0 + 0.5 * k_recomb_H_Cl * C_Cl_0 * M_0)
			/ molMass[3] };

	const scalar A45 { (-k_diss_HCl * M_0) / molMass[4] };

	const scalar B4 { C_H_0 / (timeStep * rho_0) };

	const scalar A51 { (-0.5 * k_prop_H_Cl2 * C_H_0) / molMass[0] };
	const scalar A52 { (-0.5 * k_prop_Cl_H2 * C_H2_0
			- 0.5 * k_recomb_H_Cl * C_H_0 * M_0) / molMass[1] };
	const scalar A53 { (-0.5 * k_prop_Cl_H2 * C_Cl_0) / molMass[2] };
	const scalar A54 { (-0.5 * k_prop_H_Cl2 * C_Cl2_0
			- 0.5 * k_recomb_H_Cl * C_Cl_0 * M_0) / molMass[3] };

	const scalar A55 { (1 / timeStep + k_diss_HCl * M_0) / molMass[4] };

	const scalar B5 { C_HCl_0 / (timeStep * rho_0) };

	matrixDiagonale[0] = A11;
	matrixDiagonale[1] = A22;
	matrixDiagonale[2] = A33;
	matrixDiagonale[3] = A44;
	matrixDiagonale[4] = A55;

	matrixLeftTriangle[1][0].first = A21;

	matrixLeftTriangle[2][0].first = A32;

	matrixLeftTriangle[3][0].first = A41;
	matrixLeftTriangle[3][1].first = A42;
	matrixLeftTriangle[3][2].first = A43;

	matrixLeftTriangle[4][0].first = A51;
	matrixLeftTriangle[4][1].first = A52;
	matrixLeftTriangle[4][2].first = A53;
	matrixLeftTriangle[4][3].first = A54;

	matrixRightTriangle[0][0].first = A12;
	matrixRightTriangle[0][1].first = A14;

	matrixRightTriangle[1][0].first = A23;
	matrixRightTriangle[1][1].first = A24;
	matrixRightTriangle[1][2].first = A25;

	matrixRightTriangle[2][0].first = A34;

	matrixRightTriangle[3][0].first = A45;

	matrixFreeTerm[0] = B1;
	matrixFreeTerm[1] = B2;
	matrixFreeTerm[2] = B3;
	matrixFreeTerm[3] = B4;
	matrixFreeTerm[4] = B5;
}

std::array<schemi::scalar, 5> schemi::chemicalKineticsH2Cl2Combustion::cellReactionMatrix::solve(
		const std::array<scalar, 5> & oldField,
		const std::size_t maxIterationNumber) const noexcept
{
	std::valarray<scalar> oldIteration { oldField[0], oldField[1], oldField[2],
			oldField[3], oldField[4] };

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
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3], newIteration[4]};
		}
		else if (nIterations >= maxIterationNumber)
		{
			std::clog
					<< "Gauss-Seidel algorithm for H2 + Cl2 combustion did not converged. Difference is: "
					<< diff << std::endl;

			normalize(newIteration);
			return
			{	newIteration[0], newIteration[1], newIteration[2], newIteration[3], newIteration[4]};
		}
		else
			oldIteration = newIteration;
	}
}

std::vector<schemi::chemicalKineticsH2Cl2Combustion::cellReactionMatrix> schemi::chemicalKineticsH2Cl2Combustion::velocityCalculation(
		const scalar timestep,
		const homogeneousPhase<cubicCell> & phase) const noexcept
{
	const auto size = phase.pressure.meshRef().cellsSize();

	std::vector<cellReactionMatrix> concentrationVelocityMatrix(size);

	for (std::size_t i = 0; i < size; ++i)
	{
		const scalar T = phase.temperature.ref()[i];

		const scalar k_Cl2_diss = A_Cl2_diss * std::pow(T, n_Cl2_diss)
				* std::exp(-E_Cl2_diss / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_Cl_recomb = A_Cl_recomb * std::pow(T, n_Cl_recomb)
				* std::exp(
						-E_Cl_recomb / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_H2_diss = A_H2_diss * std::pow(T, n_H2_diss)
				* std::exp(-E_H2_diss / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_H_recomb = A_H_recomb * std::pow(T, n_H_recomb)
				* std::exp(-E_H_recomb / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_Cl_H2_prop = A_Cl_H2_prop * std::pow(T, n_Cl_H2_prop)
				* std::exp(
						-E_Cl_H2_prop / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_H_Cl2_prop = A_H_Cl2_prop * std::pow(T, n_H_Cl2_prop)
				* std::exp(
						-E_H_Cl2_prop / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_H_Cl_recomb = A_H_Cl_recomb * std::pow(T, n_H_Cl_recomb)
				* std::exp(
						-E_H_Cl_recomb / (phase.phaseThermodynamics->Rv() * T));

		const scalar k_HCl_diss = A_HCl_diss * std::pow(T, n_HCl_diss)
				* std::exp(-E_HCl_diss / (phase.phaseThermodynamics->Rv() * T));

		const scalar M = phase.concentration.v[0].ref()[i];
		const scalar Cl2 = phase.concentration.v[1].ref()[i];
		const scalar Cl = phase.concentration.v[2].ref()[i];
		const scalar H2 = phase.concentration.v[3].ref()[i];
		const scalar H = phase.concentration.v[4].ref()[i];
		const scalar HCl = phase.concentration.v[5].ref()[i];

		const scalar rho = phase.density[0].ref()[i];

		concentrationVelocityMatrix[i] = cellReactionMatrix(timestep,
				k_Cl2_diss, k_Cl_recomb, k_H2_diss, k_H_recomb, k_Cl_H2_prop,
				k_H_Cl2_prop, k_H_Cl_recomb, k_HCl_diss, Cl2, Cl, H2, H, HCl, M,
				rho, phase.phaseThermodynamics->Mv());
	}

	return concentrationVelocityMatrix;
}

void schemi::chemicalKineticsH2Cl2Combustion::timeStepIntegration(
		homogeneousPhase<cubicCell> & phaseN) const noexcept
{
	auto & mesh_ = phaseN.pressure.meshRef();

	const scalar timeStep = mesh_.timestep();

	volumeField<scalar> deltaU(mesh_, 0);

	std::valarray<scalar> sourceTimeStep(veryBig, mesh_.cellsSize());

	const auto vC = velocityCalculation(timeStep, phaseN);

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const std::array<scalar, 5> oldMassFraction {
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
						/ phaseN.density[0].ref()[i],
				phaseN.concentration.v[5].ref()[i]
						* phaseN.phaseThermodynamics->Mv()[4]
						/ phaseN.density[0].ref()[i] };

		auto [newCl2, newCl, newH2, newH, newHCl] = vC[i].solve(oldMassFraction,
				maxIterationNumber);

		newCl2 = std::max(0., newCl2);
		newCl = std::max(0., newCl);
		newH2 = std::max(0., newH2);
		newH = std::max(0., newH);
		newHCl = std::max(0., newHCl);

		const auto sumFrac = newCl2 + newCl + newH2 + newH + newHCl;
		newCl2 /= sumFrac;
		newCl /= sumFrac;
		newH2 /= sumFrac;
		newH /= sumFrac;
		newHCl /= sumFrac;

		const auto newCCl2 = newCl2 * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[0];
		const auto newCCl = newCl * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[1];
		const auto newCH2 = newH2 * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[2];
		const auto newCH = newH * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[3];
		const auto newCHCl = newHCl * phaseN.density[0].ref()[i]
				/ phaseN.phaseThermodynamics->Mv()[4];

		{
			const scalar deltaC_HCl = newCHCl
					- phaseN.concentration.v[5].ref_r()[i];

			const auto & thermo = *phaseN.phaseThermodynamics;

			const auto deltaCv = 2 * thermo.Cvv()[4] - thermo.Cvv()[0]
					- thermo.Cvv()[2];

			deltaU.ref_r()[i] = -(ΔU_298
					+ deltaCv * (phaseN.temperature.ref()[i] - 298.15))
					* deltaC_HCl;
		}

		phaseN.concentration.v[1].ref_r()[i] = newCCl2;
		phaseN.concentration.v[2].ref_r()[i] = newCCl;
		phaseN.concentration.v[3].ref_r()[i] = newCH2;
		phaseN.concentration.v[4].ref_r()[i] = newCH;
		phaseN.concentration.v[5].ref_r()[i] = newCHCl;

		phaseN.concentration.v[0].ref_r()[i] = newCCl2 + newCCl + newCH2 + newCH
				+ newCHCl;
	}

	phaseN.internalEnergy.ref_r() += deltaU.ref();

	for (std::size_t k = 1; k < phaseN.density.size(); ++k)
		phaseN.density[k].ref_r() = phaseN.concentration.v[k].ref()
				* phaseN.phaseThermodynamics->Mv()[k - 1];
}

schemi::chemicalKineticsH2Cl2Combustion::chemicalKineticsH2Cl2Combustion(
		const homogeneousPhase<cubicCell> & phaseIn, const std::size_t itNumber) :
		abstractChemicalKinetics(true, itNumber)
{
	if (phaseIn.concentration.v.size() < 6)
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

	chem >> skipBuffer >> A_Cl2_diss;
	chem >> skipBuffer >> n_Cl2_diss;
	chem >> skipBuffer >> E_Cl2_diss;

	chem >> skipBuffer >> A_Cl_recomb;
	chem >> skipBuffer >> n_Cl_recomb;
	chem >> skipBuffer >> E_Cl_recomb;

	chem >> skipBuffer >> A_H2_diss;
	chem >> skipBuffer >> n_H2_diss;
	chem >> skipBuffer >> E_H2_diss;

	chem >> skipBuffer >> A_H_recomb;
	chem >> skipBuffer >> n_H_recomb;
	chem >> skipBuffer >> E_H_recomb;

	chem >> skipBuffer >> A_Cl_H2_prop;
	chem >> skipBuffer >> n_Cl_H2_prop;
	chem >> skipBuffer >> E_Cl_H2_prop;

	chem >> skipBuffer >> A_H_Cl2_prop;
	chem >> skipBuffer >> n_H_Cl2_prop;
	chem >> skipBuffer >> E_H_Cl2_prop;

	chem >> skipBuffer >> A_H_Cl_recomb;
	chem >> skipBuffer >> n_H_Cl_recomb;
	chem >> skipBuffer >> E_H_Cl_recomb;

	chem >> skipBuffer >> A_HCl_diss;
	chem >> skipBuffer >> n_HCl_diss;
	chem >> skipBuffer >> E_HCl_diss;

	chem >> skipBuffer >> ΔН_298;

	ΔU_298 = ΔН_298 - Δn * phaseIn.phaseThermodynamics->Rv() * 298.15;
}

void schemi::chemicalKineticsH2Cl2Combustion::solveChemicalKinetics(
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
