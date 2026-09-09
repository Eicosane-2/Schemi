/*
 * chemicalReactionsSystem.hpp
 *
 *  Created on: 2026/06/30
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALREACTIONSSYSTEM_HPP_
#define CHEMICALREACTIONSSYSTEM_HPP_

#include <map>
#include <vector>

#include "concentrationsPack.hpp"
#include "homogeneousPhase.hpp"
#include "irreversibleReaction.hpp"
#include "reactionMatrixCell.hpp"

namespace schemi
{
namespace chemicalKinetics
{

class chemicalReactionsSystem
{
	constexpr static scalar massFracTolerance { 1E-3 };

	bool chemicalReactions { false };

	const abstractMixtureThermodynamics & therm;

	std::vector<irreversibleReaction> reactionsParameters { };
	std::map<std::string, std::size_t> reactingComponentsMatching { };
	std::vector<std::size_t> reactingComponentsIndexes { };

	std::vector<reactionMatrixCell> matrixPrototype { };

	template<typename typeOfEntity>
	std::pair<field<std::array<std::valarray<scalar>, 2>, typeOfEntity>,
			field<scalar, typeOfEntity>> generateMatrix(
			const concentrationsPack<typeOfEntity> & c,
			const field<scalar, typeOfEntity> & T,
			const field<scalar, typeOfEntity> & rho,
			const scalar timestep) const
	{
		field<std::array<std::valarray<scalar>, 2>, typeOfEntity> Ab(
				c.v[0].meshRef(),
				{ std::valarray<scalar>(0., matrixPrototype.size()),
						std::valarray<scalar>(0.,
								std::sqrt(matrixPrototype.size())) });
		field<scalar, typeOfEntity> reactMassFracSum(c.v[0].meshRef(), 0);

		const auto rDeltat = 1.0 / timestep;

		for (std::size_t i = 0; i < Ab.size(); ++i)
		{
			scalar massFracSum { 0 };
			for (std::size_t j = 0; j < reactingComponentsIndexes.size(); ++j)
			{
				const auto k = reactingComponentsIndexes[j];

				massFracSum += c.v[k + 1].cval()[i] * therm.Mv()[k]
						/ rho.cval()[i];
			}

			auto& [A, b] = Ab.val()[i];
			for (std::size_t j = 0; j < matrixPrototype.size(); ++j)
			{
				const auto & cellData = matrixPrototype[j];
				if (!cellData.nullCell)
				{
					if (cellData.diagonalCell)
					{
						A[j] = rDeltat / cellData.molMass.first;

						b[j / reactingComponentsMatching.size()] = c.v[std::get<
								0>(cellData.indexes)].cval()[i] / rho.cval()[i]
								* rDeltat;

						for (std::size_t k = 0; k < cellData.reactParams.size();
								++k)
						{
							const auto k_r =
									std::get<0>(cellData.reactParams[k])
											* std::pow(T.cval()[i],
													std::get<1>(
															cellData.reactParams[k]))
											* std::exp(
													-std::get<2>(
															cellData.reactParams[k])
															/ (therm.Rv()
																	* T.cval()[i]));

							auto reactVel = k_r * cellData.coeffs[k]
									* cellData.reactWeight[k];

							for (const auto & rgt : cellData.comp[k])
								reactVel *= std::pow(c.v[rgt.first].cval()[i],
										rgt.second);

							A[j] += reactVel / cellData.molMass.second[k];
						}
					}
					else
					{
						for (std::size_t k = 0; k < cellData.reactParams.size();
								++k)
						{
							const auto k_r =
									std::get<0>(cellData.reactParams[k])
											* std::pow(T.cval()[i],
													std::get<1>(
															cellData.reactParams[k]))
											* std::exp(
													-std::get<2>(
															cellData.reactParams[k])
															/ (therm.Rv()
																	* T.cval()[i]));

							auto reactVel = k_r * cellData.coeffs[k]
									* cellData.reactWeight[k];

							for (const auto & rgt : cellData.comp[k])
								reactVel *= std::pow(c.v[rgt.first].cval()[i],
										rgt.second);

							A[j] += reactVel / cellData.molMass.second[k];
						}
					}
				}
			}

			reactMassFracSum.val()[i] = massFracSum;
		}

		return
		{	Ab, reactMassFracSum};
	}

	template<typename typeOfEntity>
	void timeIntegration(concentrationsPack<typeOfEntity> & concentrations,
			field<scalar, typeOfEntity> & U,
			const field<scalar, typeOfEntity> & T,
			const field<scalar, typeOfEntity> & rho,
			const scalar timestep) const
	{
		auto [AbField, reactMassFracSumOld] = generateMatrix<typeOfEntity>(
				concentrations, T, rho, timestep);

		for (std::size_t i = 0; i < rho.size(); ++i)
		{
			auto& [A, b] = AbField.val()[i];

			auto massFracs = GaussElemination(A, b);

			renormalization(reactMassFracSumOld.cval()[i], massFracs);

			scalar deltaH { 0 };
			for (std::size_t j = 0; j < massFracs.size(); ++j)
			{
				const std::size_t k = reactingComponentsIndexes[j];

				const auto concNew = massFracs[j] * rho.cval()[i]
						/ therm.Mv()[k];

				deltaH += (concNew - concentrations.v[k + 1].cval()[i])
						* therm.dHfv()[k];

				concentrations.v[k + 1].val()[i] = concNew;
			}

			U.val()[i] -= deltaH;
		}
	}

	std::valarray<scalar> GaussElemination(std::valarray<scalar> & A,
			std::valarray<scalar> & b) const;

	void renormalization(const scalar sumMassFracOld,
			std::valarray<scalar> & massFractions) const;
public:
	chemicalReactionsSystem(const abstractMixtureThermodynamics & thermIn);
	template<typename typeOfEnity>
	void solve(homogeneousPhase<typeOfEnity> & phase,
			const MPIHandler & parall) const
	{
		if (chemicalReactions)
		{
			const auto oldConcentrations = phase.concentration;
			const auto oldInternalEnergy = phase.internalEnergy;

			timeIntegration(phase.concentration, phase.internalEnergy,
					phase.temperature, phase.density[0],
					phase.temperature.meshRef().timestep());

			phase.temperature.val() = phase.phaseThermodynamics->TFromUv(
					phase.concentration.p, phase.internalEnergy.cval());

			timeIntegration(phase.concentration, phase.internalEnergy,
					phase.temperature, phase.density[0],
					phase.temperature.meshRef().timestep());

			for (std::size_t i = 0; i < reactingComponentsIndexes.size(); ++i)
			{
				const auto k = reactingComponentsIndexes[i] + 1;

				phase.concentration.v[k].val() = ((phase.concentration.v[k]
						+ oldConcentrations.v[k]) / 2).cval();

				phase.density[k].val() = phase.concentration.v[k].cval()
						* therm.Mv()[k - 1];
			}

			phase.concentration.v[0].val() = 0;
			for (std::size_t k = 1; k < phase.concentration.v.size(); ++k)
				phase.concentration.v[0].val() +=
						phase.concentration.v[k].cval();

			phase.internalEnergy.val() = ((phase.internalEnergy
					+ oldInternalEnergy) / 2).cval();

			phase.pressure.val() = phase.phaseThermodynamics->pFromUv(
					phase.concentration.p, phase.internalEnergy.cval());
			phase.temperature.val() = phase.phaseThermodynamics->TFromUv(
					phase.concentration.p, phase.internalEnergy.cval());

			if (phase.temperature.cval().min() < 0.)
				throw exception("Negative temperature after chemical reaction.",
						errors::negativeTemperatureError);

			{
				const auto v2 = phase.velocity & phase.velocity;

				phase.totalEnergy.val() = (phase.internalEnergy
						+ phase.density[0] * v2 * 0.5 + phase.rhokTurb).cval();
			}

			phase.HelmholtzEnergy.val() = phase.phaseThermodynamics->Fv(
					phase.concentration.p, phase.temperature.cval());

			phase.entropy.val() = phase.phaseThermodynamics->Sv(
					phase.concentration.p, phase.temperature.cval());

			parall.correctBoundaryValues(phase);
		}
	}
};

}  // namespace chemicalKinetics
}  // namespace schemi

#endif /* CHEMICALREACTIONSSYSTEM_HPP_ */
