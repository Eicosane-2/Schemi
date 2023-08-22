/*
 * hardSpheresTransportModel.hpp
 *
 *  Created on: 2023/08/07
 *      Author: Maxim Boldyrev
 */

#ifndef HARDSPHERESTRANSPORTMODEL_HPP_
#define HARDSPHERESTRANSPORTMODEL_HPP_

#include <valarray>

#include "abstractMixtureThermodynamics.hpp"
#include "abstractTransportModel.hpp"
#include "scalar.hpp"
#include "volumeField.hpp"
#include "concentrationsPack.hpp"

namespace schemi
{
class hardSpheresTransportModel: public abstractTransportModel
{
	std::valarray<scalar> diameters;

	scalar muHardSph(const scalar M, const scalar T,
			const scalar sigma) const noexcept;

	scalar DHardSph(const scalar M, const scalar T, const scalar p,
			const scalar sigma) const noexcept;

	scalar kappaHardSph(const scalar M, const scalar T,
			const scalar sigma) const noexcept;
public:
	hardSpheresTransportModel(const scalar mu_in, const scalar D_in,
			const scalar kappa_in,
			const std::valarray<scalar> & diameters_in) noexcept;

	volumeField<scalar> calculateMu(const std::valarray<scalar> & M,
			const volumeField<scalar> & temperature,
			const concentrationsPack<cubicCell> & concentrations) const noexcept
					override;

	volumeField<std::valarray<std::valarray<scalar>>> calculateDm(
			const std::valarray<scalar> & M,
			const volumeField<scalar> & temperature,
			const volumeField<scalar> & pressure,
			const concentrationsPack<cubicCell> & concentrations) const noexcept
					override;

	volumeField<std::valarray<scalar>> calculateDs(
			const std::valarray<scalar> & M,
			const volumeField<scalar> & temperature,
			const volumeField<scalar> & pressure,
			const concentrationsPack<cubicCell> & concentrations) const noexcept
					override;

	volumeField<scalar> calculateKappa(const std::valarray<scalar> & M,
			const volumeField<scalar> & temperature,
			const concentrationsPack<cubicCell> & concentrations) const noexcept
					override;

	surfaceField<scalar> calculateMu(const std::valarray<scalar> & M,
			const surfaceField<scalar> & temperature,
			const concentrationsPack<quadraticSurface> & concentrations) const noexcept
					override;

	surfaceField<std::valarray<std::valarray<scalar>>> calculateDm(
			const std::valarray<scalar> & M,
			const surfaceField<scalar> & temperature,
			const surfaceField<scalar> & pressure,
			const concentrationsPack<quadraticSurface> & concentrations) const noexcept
					override;

	surfaceField<std::valarray<scalar>> calculateDs(
			const std::valarray<scalar> & M,
			const surfaceField<scalar> & temperature,
			const surfaceField<scalar> & pressure,
			const concentrationsPack<quadraticSurface> & concentrations) const noexcept
					override;

	surfaceField<scalar> calculateKappa(const std::valarray<scalar> & M,
			const surfaceField<scalar> & temperature,
			const concentrationsPack<quadraticSurface> & concentrations) const noexcept
					override;
};
}  // namespace schemi

#endif /* HARDSPHERESTRANSPORTMODEL_HPP_ */
