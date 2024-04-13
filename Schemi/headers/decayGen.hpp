/*
 * decayGen.hpp
 *
 *  Created on: 2020/03/12
 *      Author: Maxim Boldyrev
 *
 *      Class for decaying model turbulence generation.
 */

#ifndef DECAYGEN_HPP_
#define DECAYGEN_HPP_

#include "abstractTurbulenceGen.hpp"

namespace schemi
{
class decayGen: public abstractTurbulenceGen
{
public:
	decayGen(const mesh & meshIn, const bool turb_in,
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
			const diffusiveFields & diffFieldsOld, const volumeField<tensor>&,
			const volumeField<vector>&, const volumeField<vector>&,
			const volumeField<vector>&, const volumeField<tensor>&,
			const volumeField<scalar>&, const volumeField<vector>&,
			const volumeField<tensor>&, const volumeField<tensor>&,
			const volumeField<vector>&, const abstractMixtureThermodynamics&,
			const volumeField<scalar>&) const noexcept override;
};
}  // namespace schemi

#endif /* DECAYGEN_HPP_ */
