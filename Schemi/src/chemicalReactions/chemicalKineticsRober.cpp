/*
 * chemicalKineticsRober.cpp
 *
 *  Created on: 2024/06/01
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsRober.hpp"

#include <iostream>
#include <fstream>
#include <numeric>

void schemi::chemicalKinetics::Rober::cellReactionMatrix::reactionMatrix::transpose() noexcept
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

schemi::chemicalKinetics::Rober::cellReactionMatrix::cellReactionMatrix() noexcept :
		solverFlag(iterativeSolver::noSolver), matrix()
{
}

schemi::chemicalKinetics::Rober::cellReactionMatrix::cellReactionMatrix(
		const iterativeSolver solverType) :
		solverFlag(solverType), matrix()
{
	if (!solverF)
	{
		switch (solverFlag)
		{
		case iterativeSolver::GaussSeidel:
			solverF = my_solveGS;
			break;
		case iterativeSolver::ConjugateGradient:
			solverF = my_solveCG;
			break;
		case iterativeSolver::JacobiConjugateGradient:
			solverF = my_solveJCG;
			break;
		case iterativeSolver::Jacobi:
			solverF = my_solveJ;
			break;
		case iterativeSolver::GaussElimination:
			solverF = my_solveGE;
			break;
		[[unlikely]] default:
			throw exception("Unknown chemical iterative solver type.",
					errors::initialisationError);
			break;
		}
	}
}

schemi::chemicalKinetics::Rober::cellReactionMatrix::cellReactionMatrix(
		const scalar timeStep, const scalar k_1, const scalar k_2,
		const scalar k_3, const scalar C_1_0, const scalar C_2_0,
		const scalar C_3_0, [[maybe_unused]] const scalar M_0,
		const scalar rho_0, const std::array<scalar, N> & molMass,
		const iterativeSolver solverType) :
		solverFlag(solverType), matrix()
{
	const scalar A11 { (1 / timeStep + k_1) / molMass[0] };
	const scalar A12 { (-0.5 * k_3 * C_3_0) / molMass[1] };
	const scalar A13 { (-0.5 * k_3 * C_2_0) / molMass[2] };
	const scalar B1 { C_1_0 / (rho_0 * timeStep) };

	const scalar A21 { (-k_1) / molMass[0] };
	const scalar A22 { (1 / timeStep + k_2 * C_2_0 + 0.5 * k_3 * C_3_0)
			/ molMass[1] };
	const scalar A23 { (0.5 * k_3 * C_2_0) / molMass[2] };
	const scalar B2 { C_2_0 / (rho_0 * timeStep) };

	const scalar A32 { (-k_3 * C_2_0) / molMass[1] };
	const scalar A33 { (1 / timeStep) / molMass[2] };
	const scalar B3 { C_3_0 / (rho_0 * timeStep) };

	std::get<0>(matrix.Diagonale) = A11;
	std::get<1>(matrix.Diagonale) = A22;
	std::get<2>(matrix.Diagonale) = A33;

	std::get<1>(matrix.LeftTriangle)[0].first = A21;

	std::get<2>(matrix.LeftTriangle)[0].first = A32;

	std::get<0>(matrix.RightTriangle)[0].first = A12;
	std::get<0>(matrix.RightTriangle)[1].first = A13;

	std::get<1>(matrix.RightTriangle)[0].first = A23;

	std::get<0>(matrix.FreeTerm) = B1;
	std::get<1>(matrix.FreeTerm) = B2;
	std::get<2>(matrix.FreeTerm) = B3;
}

auto schemi::chemicalKinetics::Rober::cellReactionMatrix::solve(
		const std::array<scalar, N> & oldField,
		const std::size_t maxIterationNumber) -> std::array<
		scalar, N>
{
	return solverF(matrix, oldField, maxIterationNumber);
}

schemi::chemicalKinetics::Rober::cellReactionMatrix schemi::chemicalKinetics::Rober::velocityCalculation(
		const scalar timestep, const scalar T,
		const std::array<scalar, N + 1> & concentrations,
		const std::array<scalar, N> & molarMasses, const scalar rho,
		const scalar R) noexcept
{
	const scalar k_1 = A_1 * std::pow(T, n_1) * std::exp(-E_1 / (R * T));

	const scalar k_2 = A_2 * std::pow(T, n_2) * std::exp(-E_2 / (R * T));

	const scalar k_3 = A_3 * std::pow(T, n_3) * std::exp(-E_3 / (R * T));

	const scalar & C_1 = std::get<1>(concentrations);
	const scalar & C_2 = std::get<2>(concentrations);
	const scalar & C_3 = std::get<3>(concentrations);
	const scalar & M = std::get<0>(concentrations);

	return cellReactionMatrix(timestep, k_1, k_2, k_3, C_1, C_2, C_3, M, rho,
			molarMasses, itSolv);
}

void schemi::chemicalKinetics::Rober::timeStepIntegration(
		homogeneousPhase<cubicCell> & phaseN)
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
			concs[k] = phaseN.concentration.v[k].cval()[i];

		for (std::size_t k = 0; k < dens.size(); ++k)
			dens[k] = phaseN.density[k + 1].cval()[i];

		const cellReactingFields oldValues(phaseN.internalEnergy.cval()[i],
				phaseN.temperature.cval()[i], concs, dens);

		cellReactingFields newValues(oldValues);

		while (true)
		{
			try
			{
				newValues = oldValues;

				for (std::size_t st = 0; st < subItNum; ++st)
				{
					const std::array<scalar, N> oldMassFraction {
							newValues.density[0] / phaseN.density[0].cval()[i],
							newValues.density[1] / phaseN.density[0].cval()[i],
							newValues.density[2] / phaseN.density[0].cval()[i] };

					const scalar sumFracOld = std::accumulate(
							oldMassFraction.cbegin(), oldMassFraction.cend(),
							0.);

					if (sumFracOld == 0.0)
						continue;

					cellReactionVel.extractMatrix(
							velocityCalculation(subTimeStep,
									phaseN.temperature.cval()[i],
									{ newValues.concentration[0],
											newValues.concentration[1],
											newValues.concentration[2],
											newValues.concentration[3] },
									{ phaseN.phaseThermodynamics->Mv()[0],
											phaseN.phaseThermodynamics->Mv()[1],
											phaseN.phaseThermodynamics->Mv()[2] },
									phaseN.density[0].cval()[i],
									phaseN.phaseThermodynamics->Rv()));

					auto reactionResult = cellReactionVel.solve(oldMassFraction,
							maxIterationNumber);

					std::transform(reactionResult.cbegin(),
							reactionResult.cend(), reactionResult.begin(),
							[](const auto & w_k) 
							{
								return std::max(static_cast<scalar>(0.), w_k);
							});

					const scalar sumFracNew = std::accumulate(
							reactionResult.cbegin(), reactionResult.cend(), 0.);

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

					for (std::size_t k = 0; k < N; ++k)
						newValues.concentration[k + 1] = reactionResult[k]
								* phaseN.density[0].cval()[i]
								/ phaseN.phaseThermodynamics->Mv()[k];

					newValues.concentration[0] = 0;
					for (std::size_t k = 1; k < newValues.concentration.size();
							++k)
						newValues.concentration[0] +=
								newValues.concentration[k];

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

		phaseN.internalEnergy.val()[i] = newValues.internalEnergy;
		phaseN.temperature.val()[i] = newValues.temperature;
		for (std::size_t k = 0; k < newValues.concentration.size(); ++k)
			phaseN.concentration.v[k].val()[i] = newValues.concentration[k];
		for (std::size_t k = 0; k < newValues.density.size(); ++k)
			phaseN.density[k + 1].val()[i] = newValues.density[k];
	}
}

schemi::chemicalKinetics::Rober::Rober(
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

	chem >> skipBuffer >> A_1;
	chem >> skipBuffer >> n_1;
	chem >> skipBuffer >> E_1;

	chem >> skipBuffer >> A_2;
	chem >> skipBuffer >> n_2;
	chem >> skipBuffer >> E_2;

	chem >> skipBuffer >> A_3;
	chem >> skipBuffer >> n_3;
	chem >> skipBuffer >> E_3;

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

	cellReactionVel = cellReactionMatrix(itSolv);

	chem >> skipBuffer >> maxIterationNumber;

	chem.close();
}

void schemi::chemicalKinetics::Rober::solveChemicalKinetics(
		homogeneousPhase<cubicCell> & phaseIn)
{
	auto phaseN1 = phaseIn;

	timeStepIntegration(phaseN1);

	timeStepIntegration(phaseN1);

	phaseN1.pressure.val() = phaseN1.phaseThermodynamics->pFromUv(
			phaseN1.concentration.p, phaseN1.internalEnergy.cval());

	phaseN1.HelmholtzEnergy.val() = phaseN1.phaseThermodynamics->Fv(
			phaseN1.concentration.p, phaseN1.temperature.cval());

	phaseN1.entropy.val() = phaseN1.phaseThermodynamics->Sv(
			phaseN1.concentration.p, phaseN1.temperature.cval());

	{
		const auto v2 = phaseN1.velocity & phaseN1.velocity;

		phaseN1.totalEnergy.val() = (phaseN1.internalEnergy
				+ phaseN1.density[0] * v2 * 0.5 + phaseN1.rhokTurb).cval();
	}

	phaseIn.average(phaseN1, *phaseIn.phaseThermodynamics, 0.5);
}
