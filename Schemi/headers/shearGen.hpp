/*
 * shearGen.hpp
 *
 *  Created on: 2020/03/15
 *      Author: Maxim Boldyrev
 *
 *      Class for k-epsilon model turbulence generation.
 */

#ifndef SHEARGEN_HPP_
#define SHEARGEN_HPP_

#include "abstractTurbulenceGen.hpp"

namespace schemi
{
class shearGen: public abstractTurbulenceGen
{
public:
	shearGen(const mesh & meshIn, const bool turb_in,
			const turbulenceModel tm_in) noexcept;

	std::tuple<std::pair<volumeField<scalar>, volumeField<scalar>>,
			std::pair<volumeField<scalar>, volumeField<scalar>>,
			std::pair<volumeField<vector>, volumeField<vector>>,
			std::pair<volumeField<scalar>, volumeField<scalar>>,
			volumeField<scalar>> calculate(scalar & sourceTimestep,
			const scalar sourceTimestepCoeff,
			const bunchOfFields<cubicCell> & cellFields,
			const diffusiveFields & diffFieldsOld,
			const volumeField<tensor> & gradV, const volumeField<vector>&,
			const volumeField<vector>&, const volumeField<vector>&,
			const volumeField<tensor>&, const volumeField<scalar>&,
			const volumeField<vector>&, const volumeField<tensor> & spherR,
			const volumeField<tensor> & devR, const volumeField<vector>&,
			const abstractMixtureThermodynamics&,
			const volumeField<scalar>&) const noexcept override;
};
}  // namespace schemi

#endif /* SHEARGEN_HPP_ */
