/*
 * BHRModel.hpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#ifndef BHRMODEL_HPP_
#define BHRMODEL_HPP_

#include "kEpsModels.hpp"
#include "diffusiveFields.hpp"

namespace schemi
{
class BHRModel: public kEpsModels
{
	scalar thetaA(const vector & a, const scalar k,
			const scalar b) const noexcept;

	scalar C0() const noexcept;
	scalar C1() const noexcept;
	scalar C2() const noexcept;
	scalar C3() const noexcept;

	scalar Ca() const noexcept;
	scalar Cb() const noexcept;

	scalar CaMax() const noexcept;

	scalar CMSA() const noexcept;
	scalar CMSM() const noexcept;

public:
	BHRModel(const mesh & meshIn, const bool turb_in);

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
			const volumeField<vector> & gradb,
			const volumeField<tensor> & spherR,
			const volumeField<tensor> & devR,
			const volumeField<vector> & gradMav_n,
			const abstractMixtureThermodynamics & mixture,
			const volumeField<scalar> & nu_t) const noexcept override;

	std::valarray<scalar> rhoepsilon(const bunchOfFields<cubicCell> & cf,
			const abstractMixtureThermodynamics & th) const noexcept override;
};
}  // namespace schemi

#endif /* BHRMODEL_HPP_ */
