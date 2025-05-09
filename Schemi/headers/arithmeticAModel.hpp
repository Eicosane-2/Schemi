/*
 * arithmeticAModel.hpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#ifndef ARITHMETICAMODEL_HPP_
#define ARITHMETICAMODEL_HPP_

#include "kEpsModels.hpp"
#include "diffusiveFields.hpp"
#include "homogeneousPhase.hpp"

namespace schemi
{
class arithmeticAModel: public kEpsModels
{
	scalar C0() const noexcept;
	scalar C1() const noexcept;
	scalar C2() const noexcept;
	scalar C3() const noexcept;

public:
	scalar sigmaRho() const noexcept;

	arithmeticAModel(const mesh & meshIn, const bool turb_in,
			const turbulenceModel model);

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

	volumeField<vector> calculate_a(const turbulenceModel model,
			const mesh & msh, const homogeneousPhase<cubicCell> & gasPhase,
			const volumeField<vector> & gradRho,
			const volumeField<vector> & gradP) const noexcept;
};
}  // namespace schemi

#endif /* ARITHMETICAMODEL_HPP_ */
