/*
 * abstractTransportModel.hpp
 *
 *  Created on: 2023/08/07
 *      Author: Maxim Boldyrev
 */

#ifndef ABSTRACTTRANSPORTMODEL_HPP_
#define ABSTRACTTRANSPORTMODEL_HPP_

#include <memory>

#include "abstractMixtureThermodynamics.hpp"
#include "concentrationsPack.hpp"
#include "scalar.hpp"
#include "volumeField.hpp"
#include "surfaceField.hpp"
#include "transportModelEnum.hpp"

namespace schemi
{
class abstractTransportModel
{
	const scalar mu, D, kappa;
protected:
	scalar muConst() const noexcept;

	scalar DConst() const noexcept;

	scalar kappaConst() const noexcept;
public:
	abstractTransportModel(const scalar mu_in, const scalar D_in,
			const scalar kappa_in) noexcept;

	virtual ~abstractTransportModel() noexcept;

	static std::unique_ptr<abstractTransportModel> createTransportModel(
			const std::vector<std::vector<std::string>> & matrixOfSubstancesConditions,
			const scalar constNu, const scalar constD, const scalar constKappa,
			const transportModel model) noexcept;

	virtual volumeField<scalar> calculateMu(const std::valarray<scalar>&,
			const volumeField<scalar>&,
			const concentrationsPack<cubicCell>&) const noexcept;

	virtual volumeField<std::valarray<std::valarray<scalar>>> calculateDm(
			const std::valarray<scalar>&, const volumeField<scalar>&,
			const volumeField<scalar>&,
			const concentrationsPack<cubicCell>&) const noexcept;

	virtual volumeField<std::valarray<scalar>> calculateDs(
			const std::valarray<scalar>&, const volumeField<scalar>&,
			const volumeField<scalar>&,
			const concentrationsPack<cubicCell>&) const noexcept;

	virtual volumeField<scalar> calculateKappa(const std::valarray<scalar>&,
			const volumeField<scalar>&,
			const concentrationsPack<cubicCell>&) const noexcept;

	virtual surfaceField<scalar> calculateMu(const std::valarray<scalar>&,
			const surfaceField<scalar>&,
			const concentrationsPack<quadraticSurface>&) const noexcept;

	virtual surfaceField<std::valarray<std::valarray<scalar>>> calculateDm(
			const std::valarray<scalar>&, const surfaceField<scalar>&,
			const surfaceField<scalar>&,
			const concentrationsPack<quadraticSurface>&) const noexcept;

	virtual surfaceField<std::valarray<scalar>> calculateDs(
			const std::valarray<scalar>&, const surfaceField<scalar>&,
			const surfaceField<scalar>&,
			const concentrationsPack<quadraticSurface>&) const noexcept;

	virtual surfaceField<scalar> calculateKappa(const std::valarray<scalar>&,
			const surfaceField<scalar>&,
			const concentrationsPack<quadraticSurface>&) const noexcept;
};
}  // namespace schemi

#endif /* ABSTRACTTRANSPORTMODEL_HPP_ */
