/*
 * decayModel.hpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#ifndef DECAYMODEL_HPP_
#define DECAYMODEL_HPP_

#include "kEpsModels.hpp"
#include "diffusiveFields.hpp"
#include "MPIHandler.hpp"

namespace schemi
{
class decayModel: public kEpsModels
{
	scalar C2() const noexcept;

	scalar Ca() const noexcept;
	scalar Cb() const noexcept;

public:
	decayModel(const mesh & meshIn, const bool turb_in);

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
};
}  // namespace schemi

#endif /* DECAYMODEL_HPP_ */
