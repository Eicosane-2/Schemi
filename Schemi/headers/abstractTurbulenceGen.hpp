/*
 * abstractTurbulenceGen.hpp
 *
 *  Created on: 2020/03/12
 *      Author: Maxim Boldyrev
 *
 *      Interface class for turbulence generation model.
 */

#ifndef ABSTRACTTURBULENCEGEN_HPP_
#define ABSTRACTTURBULENCEGEN_HPP_

#include <tuple>
#include <memory>

#include "turbulenceModelEnum.hpp"
#include "bunchOfFields.hpp"
#include "volumeField.hpp"
#include "abstractMixtureThermodynamics.hpp"
#include "abstractTurbulentParameters.hpp"
#include "cubicCell.hpp"
#include "diffusiveFields.hpp"
#include "scalar.hpp"
#include "tensor.hpp"
#include "vector.hpp"

namespace schemi
{
class abstractTurbulenceGen
{
public:
	virtual ~abstractTurbulenceGen() noexcept =0;

	abstractTurbulenceGen(const bool turb_in,
			const turbulenceModelEnum tm_in) noexcept;

	virtual std::tuple<volumeField<scalar>, volumeField<scalar>,
			volumeField<vector>, volumeField<scalar>> calculate(
			scalar& /*sourceTimestep*/, const scalar /*sourceTimestepCoeff*/,
			const bunchOfFields<cubicCell>& /*cellFields*/,
			const diffusiveFields& /*diffFieldsOld*/,
			const volumeField<tensor>& /*gradV*/,
			const volumeField<vector>& /*divDevPhysVisc*/,
			const volumeField<vector>& /*gradP*/,
			const volumeField<vector>& /*gradRho*/,
			const volumeField<tensor>& /*grada*/,
			const volumeField<scalar>& /*diva*/,
			const volumeField<vector>& /*gradb*/,
			const volumeField<tensor>& /*spherR*/,
			const volumeField<tensor>& /*devR*/,
			const volumeField<vector>& /*gradMavNormalized*/,
			const abstractMixtureThermodynamics& /*mixture*/,
			const volumeField<scalar>& /*nu_t*/) const noexcept =0;

	const bool turbulence;
	const turbulenceModelEnum model;

	std::unique_ptr<abstractTurbulentParameters> turbPar = nullptr;
};
}  // namespace schemi

#endif /* ABSTRACTTURBULENCEGEN_HPP_ */
