/*
 * BHRGen.hpp
 *
 *  Created on: 2020/03/12
 *      Author: Maxim Boldyrev
 *
 *      Class for k-epsilon-a-b model turbulence generation.
 */

#ifndef BHRGEN_HPP_
#define BHRGEN_HPP_

#include "abstractTurbulenceGen.hpp"

namespace schemi
{
class BHRGen: public abstractTurbulenceGen
{
public:
	BHRGen(const mesh & meshIn, const bool turb_in,
			const turbulenceModel tm_in) noexcept;

	std::tuple<
			std::pair<schemi::volumeField<schemi::scalar>,
					schemi::volumeField<schemi::scalar>>,
			std::pair<schemi::volumeField<schemi::scalar>,
					schemi::volumeField<schemi::scalar>>,
			std::pair<schemi::volumeField<schemi::vector>,
					schemi::volumeField<schemi::vector>>,
			std::pair<schemi::volumeField<schemi::scalar>,
					schemi::volumeField<schemi::scalar>>> calculate(
			scalar & sourceTimestep, const scalar sourceTimestepCoeff,
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

#endif /* BHRGEN_HPP_ */
