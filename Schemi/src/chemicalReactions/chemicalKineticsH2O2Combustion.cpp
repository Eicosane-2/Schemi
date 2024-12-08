/*
 * chemicalKineticsH2O2Combustion.cpp
 *
 *  Created on: 2024/02/11
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsH2O2Combustion.hpp"

#include <iostream>
#include <fstream>
#include <numeric>

void schemi::chemicalKinetics::H2O2Combustion::cellReactionMatrix::reactionMatrix::transpose() noexcept
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

schemi::chemicalKinetics::H2O2Combustion::cellReactionMatrix::cellReactionMatrix() noexcept :
		solverFlag(iterativeSolver::noSolver), matrix()
{
}

schemi::chemicalKinetics::H2O2Combustion::cellReactionMatrix::cellReactionMatrix(
		const scalar timeStep, const scalar k_R1, const scalar k_R2,
		const scalar k_R3, const scalar k_R4, const scalar k_R5,
		const scalar k_R6, const scalar k_R7, const scalar k_R8,
		const scalar k_R9, const scalar k_R10, const scalar k_R11,
		const scalar k_R12, const scalar C_O2_0, const scalar C_O_0,
		const scalar C_H2_0, const scalar C_H_0, const scalar C_OH_0,
		const scalar C_HO2_0, const scalar C_H2O_0, const scalar M_0,
		const scalar rho_0, const std::array<scalar, N> & molMass,
		const iterativeSolver solverType) :
		solverFlag(solverType), matrix()
{
	const scalar A11 { (1 / timeStep + 0.5 * (k_R1 + k_R2) * C_H2_0 + k_R4 * M_0
			+ 0.5 * k_R6 * C_H_0 + 0.5 * k_R8 * C_H_0 * M_0) / molMass[0] };
	const scalar A12 { (-k_R12 * C_O_0 * M_0) / molMass[1] };
	const scalar A13 { (0.5 * (k_R1 + k_R2) * C_O2_0) / molMass[2] };
	const scalar A14 { (0.5 * k_R6 * C_O2_0 + 0.5 * k_R8 * C_O2_0 * M_0)
			/ molMass[3] };
	const scalar B1 { C_O2_0 / (timeStep * rho_0) };

	const scalar A21 { (-2 * k_R4 * M_0 - 0.5 * k_R6 * C_H_0) / molMass[0] };
	const scalar A22 { (1 / timeStep + 0.5 * k_R5 * C_H2_0
			+ 2 * k_R12 * C_O_0 * M_0) / molMass[1] };
	const scalar A23 { (0.5 * k_R5 * C_O_0) / molMass[2] };
	const scalar A24 { (-0.5 * k_R6 * C_O2_0) / molMass[3] };
	const scalar B2 { C_O_0 / (timeStep * rho_0) };

	const scalar A31 { (0.5 * (k_R1 + k_R2) * C_H2_0) / molMass[0] };
	const scalar A32 { (0.5 * k_R5 * C_H2_0) / molMass[1] };
	const scalar A33 { (1 / timeStep + 0.5 * (k_R1 + k_R2) * C_O2_0 + k_R3 * M_0
			+ 0.5 * k_R5 * C_O_0 + 0.5 * k_R7 * C_HO2_0
			+ 0.5 * k_R9 * C_OH_0 * M_0) / molMass[2] };
	const scalar A34 { (-k_R12 * C_H_0) / molMass[3] };
	const scalar A35 { (0.5 * k_R9 * C_H2_0 * M_0) / molMass[4] };
	const scalar A36 { (0.5 * k_R7 * C_H2_0) / molMass[5] };
	const scalar B3 { C_H2_0 / (timeStep * rho_0) };

	const scalar A41 { (-0.5 * k_R2 * C_H2_0 + 0.5 * k_R6 * C_H_0
			+ 0.5 * k_R8 * C_H_0 * M_0) / molMass[0] };
	const scalar A42 { (-0.5 * k_R5 * C_H2_0) / molMass[1] };
	const scalar A43 { (-0.5 * k_R2 * C_O2_0 - 2 * k_R3 * M_0
			- 0.5 * k_R5 * C_O_0) / molMass[2] };
	const scalar A44 { (1 / timeStep + 0.5 * k_R6 * C_O2_0
			+ 0.5 * k_R8 * C_O2_0 * M_0 + k_R10 * C_H_0 * M_0
			+ 0.5 * k_R11 * C_OH_0 * M_0) / molMass[3] };
	const scalar A45 { (0.5 * k_R11 * C_H_0 * M_0) / molMass[4] };
	const scalar B4 { C_H_0 / (timeStep * rho_0) };

	const scalar A51 { (-0.5 * 2 * k_R1 * C_H2_0 - 0.5 * k_R6 * C_H_0)
			/ molMass[0] };
	const scalar A52 { (-0.5 * k_R5 * C_H2_0) / molMass[1] };
	const scalar A53 { (-0.5 * 2 * k_R1 * C_O2_0 - 0.5 * k_R5 * C_O_0
			+ 0.5 * k_R9 * C_OH_0 * M_0) / molMass[2] };
	const scalar A54 { (-0.5 * k_R6 * C_O2_0 + 0.5 * k_R11 * C_OH_0 * M_0)
			/ molMass[3] };
	const scalar A55 { (1 / timeStep + 0.5 * k_R9 * C_H2_0 * M_0
			+ 0.5 * k_R11 * C_H_0 * M_0) / molMass[4] };
	const scalar B5 { C_OH_0 / (timeStep * rho_0) };

	const scalar A61 { (-0.5 * k_R2 * C_H2_0 - 0.5 * k_R8 * C_H_0 * M_0)
			/ molMass[0] };
	const scalar A63 { (-0.5 * k_R2 * C_O2_0 + 0.5 * k_R7 * C_HO2_0)
			/ molMass[2] };
	const scalar A64 { (-0.5 * k_R8 * C_O2_0 * M_0) / molMass[3] };
	const scalar A66 { (1 / timeStep + 0.5 * k_R7 * C_H2_0) / molMass[5] };
	const scalar B6 { C_HO2_0 / (timeStep * rho_0) };

	const scalar A73 { (-0.5 * k_R7 * C_HO2_0 - 0.5 * k_R9 * C_OH_0 * M_0)
			/ molMass[2] };
	const scalar A74 { (-0.5 * k_R11 * C_OH_0 * M_0) / molMass[3] };
	const scalar A75 { (-0.5 * k_R9 * C_H2_0 * M_0 - 0.5 * k_R11 * C_H_0 * M_0)
			/ molMass[4] };
	const scalar A76 { (-0.5 * k_R7 * C_H2_0) / molMass[5] };
	const scalar A77 { (1 / timeStep) / molMass[6] };
	const scalar B7 { C_H2O_0 / (timeStep * rho_0) };

	std::get<0>(matrix.Diagonale) = A11;
	std::get<1>(matrix.Diagonale) = A22;
	std::get<2>(matrix.Diagonale) = A33;
	std::get<3>(matrix.Diagonale) = A44;
	std::get<4>(matrix.Diagonale) = A55;
	std::get<5>(matrix.Diagonale) = A66;
	std::get<6>(matrix.Diagonale) = A77;

	std::get<1>(matrix.LeftTriangle)[0].first = A21;

	std::get<2>(matrix.LeftTriangle)[0].first = A31;
	std::get<2>(matrix.LeftTriangle)[1].first = A32;

	std::get<3>(matrix.LeftTriangle)[0].first = A41;
	std::get<3>(matrix.LeftTriangle)[1].first = A42;
	std::get<3>(matrix.LeftTriangle)[2].first = A43;

	std::get<4>(matrix.LeftTriangle)[0].first = A51;
	std::get<4>(matrix.LeftTriangle)[1].first = A52;
	std::get<4>(matrix.LeftTriangle)[2].first = A53;
	std::get<4>(matrix.LeftTriangle)[3].first = A54;

	std::get<5>(matrix.LeftTriangle)[0].first = A61;
	std::get<5>(matrix.LeftTriangle)[1].first = A63;
	std::get<5>(matrix.LeftTriangle)[2].first = A64;

	std::get<6>(matrix.LeftTriangle)[0].first = A73;
	std::get<6>(matrix.LeftTriangle)[1].first = A74;
	std::get<6>(matrix.LeftTriangle)[2].first = A75;
	std::get<6>(matrix.LeftTriangle)[3].first = A76;

	std::get<0>(matrix.RightTriangle)[0].first = A12;
	std::get<0>(matrix.RightTriangle)[1].first = A13;
	std::get<0>(matrix.RightTriangle)[2].first = A14;

	std::get<1>(matrix.RightTriangle)[0].first = A23;
	std::get<1>(matrix.RightTriangle)[1].first = A24;

	std::get<2>(matrix.RightTriangle)[0].first = A34;
	std::get<2>(matrix.RightTriangle)[1].first = A35;
	std::get<2>(matrix.RightTriangle)[2].first = A36;

	std::get<3>(matrix.RightTriangle)[0].first = A45;

	std::get<0>(matrix.FreeTerm) = B1;
	std::get<1>(matrix.FreeTerm) = B2;
	std::get<2>(matrix.FreeTerm) = B3;
	std::get<3>(matrix.FreeTerm) = B4;
	std::get<4>(matrix.FreeTerm) = B5;
	std::get<5>(matrix.FreeTerm) = B6;
	std::get<6>(matrix.FreeTerm) = B7;
}

auto schemi::chemicalKinetics::H2O2Combustion::cellReactionMatrix::solve(
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
	default:
		throw exception("Unknown chemical iterative solver type.",
				errors::initialisationError);
		break;
	}
}

schemi::chemicalKinetics::H2O2Combustion::cellReactionMatrix schemi::chemicalKinetics::H2O2Combustion::velocityCalculation(
		const scalar timestep, const scalar T,
		const std::array<scalar, N + 1> & concentrations,
		const std::array<scalar, N> & molarMasses, const scalar rho,
		const scalar R) const noexcept
{
	const scalar k_R1 = A_R1 * std::pow(T, n_R1) * std::exp(-E_R1 / (R * T));
	const scalar k_R2 = A_R2 * std::pow(T, n_R2) * std::exp(-E_R2 / (R * T));
	const scalar k_R3 = A_R3 * std::pow(T, n_R3) * std::exp(-E_R3 / (R * T));
	const scalar k_R4 = A_R4 * std::pow(T, n_R4) * std::exp(-E_R4 / (R * T));
	const scalar k_R5 = A_R5 * std::pow(T, n_R5) * std::exp(-E_R5 / (R * T));
	const scalar k_R6 = A_R6 * std::pow(T, n_R6) * std::exp(-E_R6 / (R * T));
	const scalar k_R7 = A_R7 * std::pow(T, n_R7) * std::exp(-E_R7 / (R * T));
	const scalar k_R8 = A_R8 * std::pow(T, n_R8) * std::exp(-E_R8 / (R * T));
	const scalar k_R9 = A_R9 * std::pow(T, n_R9) * std::exp(-E_R9 / (R * T));
	const scalar k_R10 = A_R10 * std::pow(T, n_R10)
			* std::exp(-E_R10 / (R * T));
	const scalar k_R11 = A_R11 * std::pow(T, n_R11)
			* std::exp(-E_R11 / (R * T));
	const scalar k_R12 = A_R12 * std::pow(T, n_R12)
			* std::exp(-E_R12 / (R * T));

	const scalar & O2 = std::get<1>(concentrations);
	const scalar & O = std::get<2>(concentrations);
	const scalar & H2 = std::get<3>(concentrations);
	const scalar & H = std::get<4>(concentrations);
	const scalar & OH = std::get<5>(concentrations);
	const scalar & HO2 = std::get<6>(concentrations);
	const scalar & H2O = std::get<7>(concentrations);
	const scalar & M = std::get<0>(concentrations);

	return cellReactionMatrix(timestep, k_R1, k_R2, k_R3, k_R4, k_R5, k_R6,
			k_R7, k_R8, k_R9, k_R10, k_R11, k_R12, O2, O, H2, H, OH, HO2, H2O,
			M, rho, molarMasses, itSolv);
}

void schemi::chemicalKinetics::H2O2Combustion::timeStepIntegration(
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
							newValues.density[4] / phaseN.density[0]()[i],
							newValues.density[5] / phaseN.density[0]()[i],
							newValues.density[6] / phaseN.density[0]()[i] };

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
									newValues.concentration[5],
									newValues.concentration[6],
									newValues.concentration[7] },
							{ phaseN.phaseThermodynamics->Mv()[0],
									phaseN.phaseThermodynamics->Mv()[1],
									phaseN.phaseThermodynamics->Mv()[2],
									phaseN.phaseThermodynamics->Mv()[3],
									phaseN.phaseThermodynamics->Mv()[4],
									phaseN.phaseThermodynamics->Mv()[5],
									phaseN.phaseThermodynamics->Mv()[6] },
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
						const scalar deltaC_H2O = std::get<6>(reactionResult)
								* phaseN.density[0]()[i]
								/ phaseN.phaseThermodynamics->Mv()[6]
								- newValues.concentration[7];

						const auto & thermo = *phaseN.phaseThermodynamics;

						const auto deltaCv = 2 * thermo.Cvv()[6]
								- thermo.Cvv()[0] - 2 * thermo.Cvv()[2];

						deltaU = -(ΔU_298
								+ deltaCv * (phaseN.temperature()[i] - 298.15))
								* deltaC_H2O;
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

schemi::chemicalKinetics::H2O2Combustion::H2O2Combustion(
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
		throw std::ifstream::failure("./set/chemicalKinetics.txt not found.");

	chem >> skipBuffer >> skipBuffer;

	chem >> skipBuffer >> A_R1;
	chem >> skipBuffer >> n_R1;
	chem >> skipBuffer >> E_R1;

	chem >> skipBuffer >> A_R2;
	chem >> skipBuffer >> n_R2;
	chem >> skipBuffer >> E_R2;

	chem >> skipBuffer >> A_R3;
	chem >> skipBuffer >> n_R3;
	chem >> skipBuffer >> E_R3;

	chem >> skipBuffer >> A_R4;
	chem >> skipBuffer >> n_R4;
	chem >> skipBuffer >> E_R4;

	chem >> skipBuffer >> A_R5;
	chem >> skipBuffer >> n_R5;
	chem >> skipBuffer >> E_R5;

	chem >> skipBuffer >> A_R6;
	chem >> skipBuffer >> n_R6;
	chem >> skipBuffer >> E_R6;

	chem >> skipBuffer >> A_R7;
	chem >> skipBuffer >> n_R7;
	chem >> skipBuffer >> E_R7;

	chem >> skipBuffer >> A_R8;
	chem >> skipBuffer >> n_R8;
	chem >> skipBuffer >> E_R8;

	chem >> skipBuffer >> A_R9;
	chem >> skipBuffer >> n_R9;
	chem >> skipBuffer >> E_R9;

	chem >> skipBuffer >> A_R10;
	chem >> skipBuffer >> n_R10;
	chem >> skipBuffer >> E_R10;

	chem >> skipBuffer >> A_R11;
	chem >> skipBuffer >> n_R11;
	chem >> skipBuffer >> E_R11;

	chem >> skipBuffer >> A_R12;
	chem >> skipBuffer >> n_R12;
	chem >> skipBuffer >> E_R12;

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

	ΔU_298 = ΔН_298 - Δn * phaseIn.phaseThermodynamics->Rv() * 298.15;

	chem.close();
}

void schemi::chemicalKinetics::H2O2Combustion::solveChemicalKinetics(
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

