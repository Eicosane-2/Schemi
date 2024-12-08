/*
 * kEpsAModel.hpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#ifndef KEPSAMODEL_HPP_
#define KEPSAMODEL_HPP_

#include "kEpsModels.hpp"
#include "diffusiveFields.hpp"

namespace schemi
{
class kEpsAModel: public kEpsModels
{
	scalar thetaA(const vector & a, const scalar k,
			const scalar b) const noexcept;

	scalar C0() const noexcept;
	scalar C1() const noexcept;
	scalar C2() const noexcept;
	scalar C3() const noexcept;

	scalar Ca() const noexcept;

	scalar CaMax() const noexcept;

	scalar CMSA() const noexcept;

	scalar chi() const noexcept;

public:
	kEpsAModel(const mesh & meshIn, const bool turb_in);

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

	volumeField<scalar> calculate_b(const mesh & m,
			const std::vector<field<scalar, cubicCell>> & c,
			const std::vector<field<scalar, cubicCell>> & rho) const noexcept;
};
}  // namespace schemi

#endif /* KEPSAMODEL_HPP_ */
