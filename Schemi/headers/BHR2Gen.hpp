/*
 * BHR2Gen.hpp
 *
 *  Created on: 2024/05/01
 *      Author: Maxim Boldyrev
 *
 *      Class for k-epsilon-a-b model turbulence generation. Second version.
 */

#ifndef BHR2GEN_HPP_
#define BHR2GEN_HPP_

#include "abstractTurbulenceGen.hpp"

namespace schemi
{
class BHR2Gen: public abstractTurbulenceGen
{
public:
	BHR2Gen(const mesh & meshIn, const bool turb_in,
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
			const volumeField<vector> & gradP,
			const volumeField<vector> & gradRho,
			const volumeField<tensor> & grada, const volumeField<scalar> & diva,
			const volumeField<vector>&, const volumeField<tensor> & spherR,
			const volumeField<tensor> & devR,
			const volumeField<vector> & gradMav_n,
			const abstractMixtureThermodynamics & mixture,
			const volumeField<scalar> & nu_t) const noexcept override;
};
}  // namespace schemi

#endif /* BHR2GEN_HPP_ */
