/*
 * BHR3Model.hpp
 *
 *  Created on: 2025/03/25
 *      Author: Maxim Boldyrev
 */

#ifndef BHR3MODEL_HPP_
#define BHR3MODEL_HPP_

#include "kEpsModels.hpp"
#include "diffusiveFields.hpp"

namespace schemi
{
class BHR3Model: public kEpsModels
{
	scalar thetaA(const vector & a, const scalar k,
			const scalar b) const noexcept;

	scalar thetaB(const vector & gradMav_n, const vector & gradRho,
			const scalar pressure, const scalar k,
			const scalar eps) const noexcept;

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

	scalar Camp() const noexcept;
	scalar CMSB() const noexcept;

public:
	BHR3Model(const mesh & meshIn, const bool turb_in);

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

	std::valarray<scalar> rhoepsilon(const bunchOfFields<cubicCell> & cf,
			const abstractMixtureThermodynamics & th) const noexcept override;
};
}  // namespace schemi

#endif /* BHR3MODEL_HPP_ */
