/*
 * homogeneousPhase.hpp
 *
 *  Created on: 2022/10/09
 *      Author: Maxim Boldyrev
 */

#ifndef HOMOGENEOUSPHASE_HPP_
#define HOMOGENEOUSPHASE_HPP_

#include <memory>

#include "abstractTurbulenceModel.hpp"
#include "abstractTransportModel.hpp"
#include "bunchOfFields.hpp"
#include "transportCoefficients.hpp"

namespace schemi
{
template<typename typeOfEnity>
struct homogeneousPhase: bunchOfFields<typeOfEnity>, transportCoefficients<
		typeOfEnity>
{
	std::shared_ptr<abstractMixtureThermodynamics> phaseThermodynamics;
	std::shared_ptr<abstractTurbulenceModel> turbulence;
	std::shared_ptr<abstractTransportModel> transportModel;

	field<scalar, typeOfEnity> alphaFrac;

	homogeneousPhase(const bunchOfFields<typeOfEnity> & bunchOfFields_in,
			const transportCoefficients<typeOfEnity> & transportCoefficients_in,
			std::unique_ptr<abstractMixtureThermodynamics> & abstractMixtureThermodynamics_in,
			std::unique_ptr<abstractTurbulenceModel> & turbulence_in,
			std::unique_ptr<abstractTransportModel> & transportModel_in) noexcept :
			bunchOfFields<typeOfEnity>(bunchOfFields_in), transportCoefficients<
					typeOfEnity>(transportCoefficients_in), phaseThermodynamics(
					std::move(abstractMixtureThermodynamics_in)), turbulence(
					std::move(turbulence_in)), transportModel(
					std::move(transportModel_in)), alphaFrac(
					this->pressure.meshRef(), 1.)
	{
	}

	homogeneousPhase(const bunchOfFields<typeOfEnity> & bunchOfFields_in,
			const transportCoefficients<typeOfEnity> & transportCoefficients_in,
			std::shared_ptr<abstractMixtureThermodynamics> & abstractMixtureThermodynamics_in,
			std::shared_ptr<abstractTurbulenceModel> & turbulence_in,
			std::shared_ptr<abstractTransportModel> & transportModel_in) noexcept :
			bunchOfFields<typeOfEnity>(bunchOfFields_in), transportCoefficients<
					typeOfEnity>(transportCoefficients_in), phaseThermodynamics(
					abstractMixtureThermodynamics_in), turbulence(
					turbulence_in), transportModel(transportModel_in), alphaFrac(
					this->pressure.meshRef(), 1.)
	{
	}
};
}  // namespace schemi

#endif /* HOMOGENEOUSPHASE_HPP_ */
