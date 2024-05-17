/*
 * zeroGen.hpp
 *
 *  Created on: 2020/03/12
 *      Author: Maxim Boldyrev
 *
 *      Class for zero turbulence generation.
 */

#ifndef ZEROGEN_HPP_
#define ZEROGEN_HPP_

#include "abstractTurbulenceGen.hpp"

namespace schemi
{
class zeroGen: public abstractTurbulenceGen
{
public:
	zeroGen(const mesh & meshIn, const bool turb_in,
			const turbulenceModel tm_in) noexcept;

	std::tuple<std::pair<volumeField<scalar>, volumeField<scalar>>,
			std::pair<volumeField<scalar>, volumeField<scalar>>,
			std::pair<volumeField<vector>, volumeField<vector>>,
			std::pair<volumeField<scalar>, volumeField<scalar>>,
			volumeField<scalar>> calculate(scalar & sourceTimestep,
			const scalar, const bunchOfFields<cubicCell> & cellFields,
			const diffusiveFields&, const volumeField<tensor>&,
			const volumeField<vector>&, const volumeField<vector>&,
			const volumeField<vector>&, const volumeField<tensor>&,
			const volumeField<scalar>&, const volumeField<vector>&,
			const volumeField<tensor>&, const volumeField<tensor>&,
			const volumeField<vector>&, const abstractMixtureThermodynamics&,
			const volumeField<scalar>&) const noexcept override;
};
}  // namespace schemi

#endif /* ZEROGEN_HPP_ */
