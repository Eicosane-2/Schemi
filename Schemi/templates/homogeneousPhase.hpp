/*
 * homogeneousPhase.hpp
 *
 *  Created on: 2022/10/09
 *      Author: Maxim Boldyrev
 */

#ifndef HOMOGENEOUSPHASE_HPP_
#define HOMOGENEOUSPHASE_HPP_

#include <memory>

#include "abstractTurbulenceGen.hpp"
#include "abstractTransportModel.hpp"
#include "bunchOfFields.hpp"
#include "transportCoefficients.hpp"

namespace schemi
{
template<typename typeOfEnity1>
struct homogeneousPhase: bunchOfFields<typeOfEnity1>, transportCoefficients<
		typeOfEnity1>
{
	std::shared_ptr<abstractMixtureThermodynamics> phaseThermodynamics;
	std::shared_ptr<abstractTurbulenceGen> turbulenceSources;
	std::shared_ptr<abstractTransportModel> transportModel;

	field<scalar, typeOfEnity1> alphaFrac;

	homogeneousPhase(const bunchOfFields<typeOfEnity1> & bunchOfFields_in,
			const transportCoefficients<typeOfEnity1> & transportCoefficients_in,
			std::unique_ptr<abstractMixtureThermodynamics> & abstractMixtureThermodynamics_in,
			std::unique_ptr<abstractTurbulenceGen> & turbulenceSources_in,
			std::unique_ptr<abstractTransportModel> & transportModel_in) noexcept :
			bunchOfFields<typeOfEnity1>(bunchOfFields_in), transportCoefficients<
					typeOfEnity1>(transportCoefficients_in), phaseThermodynamics(
					std::move(abstractMixtureThermodynamics_in)), turbulenceSources(
					std::move(turbulenceSources_in)), transportModel(
					std::move(transportModel_in)), alphaFrac(
					this->pressure.meshRef(), 1.)
	{
	}

	homogeneousPhase(const bunchOfFields<typeOfEnity1> & bunchOfFields_in,
			const transportCoefficients<typeOfEnity1> & transportCoefficients_in,
			std::shared_ptr<abstractMixtureThermodynamics> & abstractMixtureThermodynamics_in,
			std::shared_ptr<abstractTurbulenceGen> & turbulenceSources_in,
			std::shared_ptr<abstractTransportModel> & transportModel_in) noexcept :
			bunchOfFields<typeOfEnity1>(bunchOfFields_in), transportCoefficients<
					typeOfEnity1>(transportCoefficients_in), phaseThermodynamics(
					abstractMixtureThermodynamics_in), turbulenceSources(
					turbulenceSources_in), transportModel(transportModel_in), alphaFrac(
					this->pressure.meshRef(), 1.)
	{
	}
};
}  // namespace schemi

#endif /* HOMOGENEOUSPHASE_HPP_ */
