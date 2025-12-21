/*
 * LennardJonesTransportModel.hpp
 *
 *  Created on: 2025/11/27
 *      Author: Maxim Boldyrev
 */

#ifndef LENNARDJONESTRANSPORTMODEL_HPP_
#define LENNARDJONESTRANSPORTMODEL_HPP_

#include <valarray>

#include "abstractMixtureThermodynamics.hpp"
#include "abstractTransportModel.hpp"
#include "scalar.hpp"
#include "volumeField.hpp"
#include "concentrationsPack.hpp"

namespace schemi
{
class LennardJonesTransportModel: public abstractTransportModel
{
	std::valarray<scalar> epsilonLJk;
	std::valarray<scalar> sigmaLJ;

	scalar OmegaCalcMuKappa(const scalar Ts) const noexcept;
	scalar OmegaCalcD(const scalar Ts) const noexcept;

	scalar muLJ(const scalar M, const scalar T, const scalar Omega,
			const scalar sigma) const noexcept;

	scalar DLJ(const scalar M, const scalar T, const scalar p,
			const scalar Omega, const scalar sigma) const noexcept;

	scalar kappaLJ(const scalar M, const scalar T, const scalar Omega,
			const scalar sigma) const noexcept;
public:
	LennardJonesTransportModel(const scalar mu_in, const scalar D_in,
			const scalar kappa_in, const std::valarray<scalar> & epsilonLJk_in,
			const std::valarray<scalar> & sigmaLJ_in) noexcept;

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

#endif /* LENNARDJONESTRANSPORTMODEL_HPP_ */
