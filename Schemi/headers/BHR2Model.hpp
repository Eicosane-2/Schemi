/*
 * BHR2Model.hpp
 *
 *  Created on: 2024/12/03
 *      Author: Maxim Boldyrev
 */

#ifndef BHR2MODEL_HPP_
#define BHR2MODEL_HPP_

#include "kEpsModels.hpp"
#include "diffusiveFields.hpp"

namespace schemi
{
class BHR2Model: public kEpsModels
{
	scalar thetaA(const vector & a, const scalar k,
			const scalar b) const noexcept;

	scalar C0() const noexcept;
	scalar C1() const noexcept;
	scalar C2() const noexcept;
	scalar C3() const noexcept;
	scalar C4() const noexcept;

	scalar Ca() const noexcept;
	scalar Cb() const noexcept;

	scalar alpha2() const noexcept;
	scalar alpha3() const noexcept;
	scalar alpha4() const noexcept;

	scalar CaMax() const noexcept;

	scalar CMSA() const noexcept;
	scalar CMSM() const noexcept;

public:
	BHR2Model(const mesh & meshIn, const bool turb_in);

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

#endif /* BHR2MODEL_HPP_ */
