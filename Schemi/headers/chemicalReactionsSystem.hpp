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

	scalar minTimeStep { 0 };

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
								reactingComponentsMatching.size()) });
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
			field<scalar, typeOfEntity> & U, field<scalar, typeOfEntity> & T,
			const field<scalar, typeOfEntity> & rho,
			const scalar timestep) const
	{
		std::size_t subSteps = 1;
		scalar subTimestep = timestep;
		while (true)
		{
			try
			{
				for (std::size_t st = 0; st < subSteps; ++st)
				{
					auto [AbField, reactMassFracSumOld] = generateMatrix<
							typeOfEntity>(concentrations, T, rho, subTimestep);

					for (std::size_t i = 0; i < rho.size(); ++i)
					{
						auto& [A, b] = AbField.val()[i];

						auto massFracs = GaussElemination(A, b);

						renormalization(reactMassFracSumOld.cval()[i],
								massFracs);

						scalar deltaCp { 0 };
						scalar deltaH { 0 };
						const auto TOld = T.cval()[i];
						for (std::size_t j = 0; j < massFracs.size(); ++j)
						{
							const std::size_t k = reactingComponentsIndexes[j];

							const auto concOld =
									concentrations.v[k + 1].cval()[i];
							const auto concNew = massFracs[j] * rho.cval()[i]
									/ therm.Mv()[k];
							const auto deltaC = concNew - concOld;

							deltaH += deltaC * therm.dHfv()[k];

							deltaCp += deltaC * therm.Cpk(concOld, TOld, k);

							concentrations.v[k + 1].val()[i] = concNew;
						}

						U.val()[i] -= (deltaH + deltaCp * (TOld - 298.15));
					}

					T.val() = therm.TFromUv(concentrations.p, U.cval());

					if (T.cval().min() < 0.)
						throw exception(
								"Negative temperature after chemical reaction.",
								errors::negativeTemperatureError);
				}
				break;
			} catch (const exception & e)
			{
				if (e.errType == errors::negativeTemperatureError
						|| e.errType == errors::positivnessError)
				{
					subSteps *= 10;
					subTimestep /= 10;

					if (subTimestep <= minTimeStep)
						throw exception(
								"Time step for chemical reactions became too small.",
								errors::systemError);

					std::cout
							<< "Chemical reaction time-step diminished. Time-step is "
							<< subTimestep << '.' << std::endl;
				}
				else
					throw e;
			}
		}
	}

	std::valarray<scalar> GaussElemination(std::valarray<scalar> & A,
			std::valarray<scalar> & b) const;

	void renormalization(const scalar sumMassFracOld,
			std::valarray<scalar> & massFractions) const;
public:
	chemicalReactionsSystem(const abstractMixtureThermodynamics & thermIn,
			const scalar minTime);
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

			{
				const auto v2 = phase.velocity & phase.velocity;

				phase.totalEnergy.val() = (phase.internalEnergy
						+ 0.5 * phase.density[0] * v2 + phase.rhokTurb).cval();
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
