/*
 * chemicalKineticsH2Cl2Combustion.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsH2Cl2Combustion.hpp"

#include <iostream>
#include <fstream>
#include <numeric>

void schemi::chemicalKinetics::H2Cl2Combustion::cellReactionMatrix::reactionMatrix::transpose() noexcept
{
	std::array<triangleList, N> LeftTriangleNew, RightTriangleNew;

	for (std::size_t i = 0; i < Diagonale.size(); ++i)
	{
		for (std::size_t j = 0; j < LeftTriangle[i].size(); ++j)
		{
			const std::size_t jAbsOld = LeftTriangle[i][j].second;
			const auto Avalue = LeftTriangle[i][j].first;

			const std::size_t iNew = jAbsOld;
			const std::size_t jNew = i;

			RightTriangleNew[iNew].emplace_back(Avalue, jNew);
		}

		for (std::size_t j = 0; j < RightTriangle[i].size(); ++j)
		{
			const std::size_t jAbsOld = RightTriangle[i][j].second;
			const auto Avalue = RightTriangle[i][j].first;

			const std::size_t iNew = jAbsOld;
			const std::size_t jNew = i;

			LeftTriangleNew[iNew].emplace_back(Avalue, jNew);
		}
	}

	LeftTriangle = LeftTriangleNew;
	RightTriangle = RightTriangleNew;
}

schemi::chemicalKinetics::H2Cl2Combustion::cellReactionMatrix::cellReactionMatrix() noexcept :
		solverFlag(iterativeSolver::noSolver), matrix()
{
}

schemi::chemicalKinetics::H2Cl2Combustion::cellReactionMatrix::cellReactionMatrix(
		const scalar timeStep, const scalar k_diss_Cl2,
		const scalar k_recomb_Cl, const scalar k_diss_H2,
		const scalar k_recomb_H, const scalar k_prop_Cl_H2,
		const scalar k_prop_H_Cl2, const scalar k_recomb_H_Cl,
		const scalar k_diss_HCl, const scalar C_Cl2_0, const scalar C_Cl_0,
		const scalar C_H2_0, const scalar C_H_0, const scalar C_HCl_0,
		const scalar M_0, const scalar rho_0,
		const std::array<scalar, N> & molMass, const iterativeSolver solverType) :
		solverFlag(solverType), matrix()
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

	std::get<0>(matrix.Diagonale) = A11;
	std::get<1>(matrix.Diagonale) = A22;
	std::get<2>(matrix.Diagonale) = A33;
	std::get<3>(matrix.Diagonale) = A44;
	std::get<4>(matrix.Diagonale) = A55;

	std::get<1>(matrix.LeftTriangle)[0].first = A21;

	std::get<2>(matrix.LeftTriangle)[0].first = A32;

	std::get<3>(matrix.LeftTriangle)[0].first = A41;
	std::get<3>(matrix.LeftTriangle)[1].first = A42;
	std::get<3>(matrix.LeftTriangle)[2].first = A43;

	std::get<4>(matrix.LeftTriangle)[0].first = A51;
	std::get<4>(matrix.LeftTriangle)[1].first = A52;
	std::get<4>(matrix.LeftTriangle)[2].first = A53;
	std::get<4>(matrix.LeftTriangle)[3].first = A54;

	std::get<0>(matrix.RightTriangle)[0].first = A12;
	std::get<0>(matrix.RightTriangle)[1].first = A14;

	std::get<1>(matrix.RightTriangle)[0].first = A23;
	std::get<1>(matrix.RightTriangle)[1].first = A24;
	std::get<1>(matrix.RightTriangle)[2].first = A25;

	std::get<2>(matrix.RightTriangle)[0].first = A34;

	std::get<3>(matrix.RightTriangle)[0].first = A45;

	std::get<0>(matrix.FreeTerm) = B1;
	std::get<1>(matrix.FreeTerm) = B2;
	std::get<2>(matrix.FreeTerm) = B3;
	std::get<3>(matrix.FreeTerm) = B4;
	std::get<4>(matrix.FreeTerm) = B5;
}

auto schemi::chemicalKinetics::H2Cl2Combustion::cellReactionMatrix::solve(
		const std::array<scalar, N> & oldField,
		const std::size_t maxIterationNumber) const -> std::array<
		scalar, N>
{
	switch (solverFlag)
	{
	case iterativeSolver::GaussSeidel:
		return solveGS<reactionMatrix, N>(matrix, oldField, maxIterationNumber);
		break;
	case iterativeSolver::ConjugateGradient:
		return solveCG<reactionMatrix, N>(matrix, oldField, maxIterationNumber);
		break;
	case iterativeSolver::JacobiConjugateGradient:
		return solveJCG<reactionMatrix, N>(matrix, oldField, maxIterationNumber);
		break;
	case iterativeSolver::Jacobi:
		return solveJ<reactionMatrix, N>(matrix, oldField, maxIterationNumber);
		break;
	case iterativeSolver::GaussElimination:
		return solveGE<reactionMatrix, N>(matrix);
		break;
	[[unlikely]] default:
		throw exception("Unknown chemical iterative solver type.",
				errors::initialisationError);
		break;
	}
}

schemi::chemicalKinetics::H2Cl2Combustion::cellReactionMatrix schemi::chemicalKinetics::H2Cl2Combustion::velocityCalculation(
		const scalar timestep, const scalar T,
		const std::array<scalar, N + 1> & concentrations,
		const std::array<scalar, N> & molarMasses, const scalar rho,
		const scalar R) const noexcept
{
	const scalar k_Cl2_diss = A_Cl2_diss * std::pow(T, n_Cl2_diss)
			* std::exp(-E_Cl2_diss / (R * T));

	const scalar k_Cl_recomb = A_Cl_recomb * std::pow(T, n_Cl_recomb)
			* std::exp(-E_Cl_recomb / (R * T));

	const scalar k_H2_diss = A_H2_diss * std::pow(T, n_H2_diss)
			* std::exp(-E_H2_diss / (R * T));

	const scalar k_H_recomb = A_H_recomb * std::pow(T, n_H_recomb)
			* std::exp(-E_H_recomb / (R * T));

	const scalar k_Cl_H2_prop = A_Cl_H2_prop * std::pow(T, n_Cl_H2_prop)
			* std::exp(-E_Cl_H2_prop / (R * T));

	const scalar k_H_Cl2_prop = A_H_Cl2_prop * std::pow(T, n_H_Cl2_prop)
			* std::exp(-E_H_Cl2_prop / (R * T));

	const scalar k_H_Cl_recomb = A_H_Cl_recomb * std::pow(T, n_H_Cl_recomb)
			* std::exp(-E_H_Cl_recomb / (R * T));

	const scalar k_HCl_diss = A_HCl_diss * std::pow(T, n_HCl_diss)
			* std::exp(-E_HCl_diss / (R * T));

	const scalar & Cl2 = std::get<1>(concentrations);
	const scalar & Cl = std::get<2>(concentrations);
	const scalar & H2 = std::get<3>(concentrations);
	const scalar & H = std::get<4>(concentrations);
	const scalar & HCl = std::get<5>(concentrations);
	const scalar & M = std::get<0>(concentrations);

	return cellReactionMatrix(timestep, k_Cl2_diss, k_Cl_recomb, k_H2_diss,
			k_H_recomb, k_Cl_H2_prop, k_H_Cl2_prop, k_H_Cl_recomb, k_HCl_diss,
			Cl2, Cl, H2, H, HCl, M, rho, molarMasses, itSolv);
}

void schemi::chemicalKinetics::H2Cl2Combustion::timeStepIntegration(
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
			concs[k] = phaseN.concentration.v[k]()[i];

		for (std::size_t k = 0; k < dens.size(); ++k)
			dens[k] = phaseN.density[k + 1]()[i];

		const cellReactingFields oldValues(phaseN.internalEnergy()[i],
				phaseN.temperature()[i], concs, dens);

		cellReactingFields newValues(oldValues);

		while (true)
		{
			try
			{
				newValues = oldValues;

				for (std::size_t st = 0; st < subItNum; ++st)
				{
					scalar deltaU { 0 };

					const std::array<scalar, N> oldMassFraction {
							newValues.density[0] / phaseN.density[0]()[i],
							newValues.density[1] / phaseN.density[0]()[i],
							newValues.density[2] / phaseN.density[0]()[i],
							newValues.density[3] / phaseN.density[0]()[i],
							newValues.density[4] / phaseN.density[0]()[i] };

					const scalar sumFracOld = std::accumulate(
							oldMassFraction.begin(), oldMassFraction.end(), 0.);

					if (sumFracOld == 0.0)
						continue;

					const auto cellReactionVel = velocityCalculation(
							subTimeStep, phaseN.temperature()[i],
							{ newValues.concentration[0],
									newValues.concentration[1],
									newValues.concentration[2],
									newValues.concentration[3],
									newValues.concentration[4],
									newValues.concentration[5] },
							{ phaseN.phaseThermodynamics->Mv()[0],
									phaseN.phaseThermodynamics->Mv()[1],
									phaseN.phaseThermodynamics->Mv()[2],
									phaseN.phaseThermodynamics->Mv()[3],
									phaseN.phaseThermodynamics->Mv()[4] },
							phaseN.density[0]()[i],
							phaseN.phaseThermodynamics->Rv());

					auto reactionResult = cellReactionVel.solve(oldMassFraction,
							maxIterationNumber);

					for (auto & w_k : reactionResult)
						w_k = std::max(static_cast<scalar>(0.), w_k);

					const scalar sumFracNew = std::accumulate(
							reactionResult.begin(), reactionResult.end(), 0.);

					if (std::abs(sumFracNew - sumFracOld) > massFracTolerance)
						throw exception(
								std::string(
										"Difference in mass fraction is too big. Delta is ")
										+ std::to_string(
												std::abs(
														sumFracNew
																- sumFracOld))
										+ '.', errors::systemError);

					for (auto & w_k : reactionResult)
					{
						w_k /= sumFracNew;
						w_k *= sumFracOld;
					}

					{
						const scalar deltaC_HCl = std::get<4>(reactionResult)
								* phaseN.density[0]()[i]
								/ phaseN.phaseThermodynamics->Mv()[4]
								- newValues.concentration[5];

						const auto & thermo = *phaseN.phaseThermodynamics;

						const auto deltaCv = 2 * thermo.Cvv()[4]
								- thermo.Cvv()[0] - thermo.Cvv()[2];

						deltaU = -(ΔU_298
								+ deltaCv * (phaseN.temperature()[i] - 298.15))
								* deltaC_HCl;
					}

					for (std::size_t k = 0; k < N; ++k)
						newValues.concentration[k + 1] = reactionResult[k]
								* phaseN.density[0]()[i]
								/ phaseN.phaseThermodynamics->Mv()[k];

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

		phaseN.internalEnergy.r()[i] = newValues.internalEnergy;
		phaseN.temperature.r()[i] = newValues.temperature;
		for (std::size_t k = 0; k < newValues.concentration.size(); ++k)
			phaseN.concentration.v[k].r()[i] = newValues.concentration[k];
		for (std::size_t k = 0; k < newValues.density.size(); ++k)
			phaseN.density[k + 1].r()[i] = newValues.density[k];
	}
}

schemi::chemicalKinetics::H2Cl2Combustion::H2Cl2Combustion(
		const homogeneousPhase<cubicCell> & phaseIn, const scalar mt) :
		abstractChemicalKinetics(true, mt), itSolv(iterativeSolver::noSolver)
{
	if (phaseIn.concentration.v.size() < N + 1)
		throw exception("Wrong number of substances.",
				errors::initialisationError);

	std::string skipBuffer;

	std::ifstream chem { "./set/chemicalKinetics.txt" };

	if (chem.is_open())
		std::cout << "./set/chemicalKinetics.txt is opened." << std::endl;
	else
		[[unlikely]]
		throw std::ifstream::failure("./set/chemicalKinetics.txt not found.");

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

	std::string solverName;

	chem >> skipBuffer >> solverName;

	try
	{
		itSolv = solvers.at(solverName);
	} catch (const std::out_of_range&)
	{
		throw exception("Unknown type of chemical iterative solver.",
				errors::initialisationError);
	}

	chem >> skipBuffer >> maxIterationNumber;

	//ΔU_298 = ΔН_298 - Δn * phaseIn.phaseThermodynamics->Rv() * 298.15;
	ΔU_298 = ΔН_298;

	chem.close();
}

void schemi::chemicalKinetics::H2Cl2Combustion::solveChemicalKinetics(
		homogeneousPhase<cubicCell> & phaseIn) const
{
	auto phaseN1 = phaseIn;

	timeStepIntegration(phaseN1);

	timeStepIntegration(phaseN1);

	phaseN1.pressure.r() = phaseN1.phaseThermodynamics->pFromUv(
			phaseN1.concentration.p, phaseN1.internalEnergy());

	phaseN1.HelmholtzEnergy.r() = phaseN1.phaseThermodynamics->Fv(
			phaseN1.concentration.p, phaseN1.temperature());

	phaseN1.entropy.r() = phaseN1.phaseThermodynamics->Sv(
			phaseN1.concentration.p, phaseN1.temperature());

	{
		const auto v2 = ampProduct(phaseN1.velocity, phaseN1.velocity);

		phaseN1.totalEnergy.r() = phaseN1.internalEnergy()
				+ phaseN1.density[0]() * v2() * 0.5 + phaseN1.rhokTurb();
	}

	phaseIn.average(phaseN1, *phaseIn.phaseThermodynamics, 0.5);
}
