/*
 * arithmeticAGen.hpp
 *
 *  Created on: 2020/03/30
 *      Author: Maxim Boldyrev
 *
 *      Class for k-epsilon-a model turbulence generation.
 */

#ifndef ARITHMETICAGEN_HPP_
#define ARITHMETICAGEN_HPP_

#include "abstractTurbulenceGen.hpp"

namespace schemi
{
class arithmeticAGen: public abstractTurbulenceGen
{
public:
	arithmeticAGen(const mesh & meshIn, const bool turb_in,
			const turbulenceModel tm_in) noexcept;

	std::tuple<std::pair<volumeField<scalar>, volumeField<scalar>>,
			std::pair<volumeField<scalar>, volumeField<scalar>>,
			std::pair<volumeField<vector>, volumeField<vector>>,
			std::pair<volumeField<scalar>, volumeField<scalar>>,
			volumeField<scalar>> calculate(scalar & sourceTimestep,
			const scalar sourceTimestepCoeff,
			const bunchOfFields<cubicCell> & cellFields,
			const diffusiveFields & diffFieldsOld,
			const volumeField<tensor> & gradV,
			const volumeField<vector> & divDevPhysVisc,
			const volumeField<vector> & gradP, const volumeField<vector>&,
			const volumeField<tensor> & grada, const volumeField<scalar> & diva,
			const volumeField<vector>&, const volumeField<tensor> & spherR,
			const volumeField<tensor> & devR, const volumeField<vector>&,
			const abstractMixtureThermodynamics&,
			const volumeField<scalar>&) const noexcept override;
};
}  // namespace schemi

#endif /* ARITHMETICAGEN_HPP_ */
