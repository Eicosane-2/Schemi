/*
 * mixtureIdeal.hpp
 *
 *  Created on: 2020/02/18
 *      Author: Maxim Boldyrev
 *
 *      Class for ideal fluids mixture thermodynamics.
 */

#ifndef MIXTUREIDEAL_HPP_
#define MIXTUREIDEAL_HPP_

#include "abstractMixtureThermodynamics.hpp"
#include "idealFluid.hpp"

namespace schemi
{
class mixtureIdeal: private idealFluid, public abstractMixtureThermodynamics
{
	const std::valarray<scalar> M, CvArr, molecMass;
public:
	mixtureIdeal() noexcept;

	mixtureIdeal(const scalar Rin, const scalar hPin,
			const std::valarray<scalar> & Min,
			const std::valarray<scalar> & Cvin) noexcept;

	scalar Rv() const noexcept override;

	const std::valarray<scalar>& Mv() const noexcept override;

	const std::valarray<scalar>& Cvv() const noexcept override;

	/**** For field ****/
	std::valarray<scalar> Cv(
			const std::vector<const std::valarray<scalar>*> & concentrations) const noexcept
					override;

	std::valarray<scalar> pFromUv(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar> & Uv) const noexcept override;

	std::valarray<scalar> UvFromp(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar> & p) const noexcept override;

	std::valarray<scalar> pcFromT(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar> & T) const noexcept override;

	std::valarray<scalar> pcFromTk(const std::valarray<scalar> & concentration,
			const std::valarray<scalar> & T, const std::size_t) const noexcept
					override;

	std::valarray<scalar> UvcFromT(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar> & T) const noexcept override;

	std::valarray<scalar> UvcFromTk(const std::valarray<scalar> & concentration,
			const std::valarray<scalar> & T,
			const std::size_t componentIndex) const noexcept override;

	std::valarray<scalar> TFromUv(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar> & Uv) const noexcept override;

	std::valarray<scalar> dpdrho(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar>&) const noexcept override;

	std::valarray<scalar> dpdUv(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar>&) const noexcept override;

	std::valarray<scalar> nonIdeality(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar>&) const noexcept override;

	std::valarray<scalar> sqSonicSpeed(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar> & density,
			const std::valarray<scalar> & Uv,
			const std::valarray<scalar> & pressure) const noexcept override;

	std::valarray<scalar> Cp(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar>&) const noexcept override;

	std::valarray<scalar> Cpk(const std::valarray<scalar>&,
			const std::valarray<scalar> & T,
			const std::size_t componentIndex) const noexcept override;

	std::valarray<scalar> hT(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar> & T) const noexcept override;

	std::valarray<scalar> hkT(const std::valarray<scalar> & concentration,
			const std::valarray<scalar> & T,
			const std::size_t componentIndex) const noexcept override;

	std::valarray<scalar> Fv(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar> & T) const noexcept override;

	std::valarray<scalar> Fvk(const std::valarray<scalar> & concentration,
			const std::valarray<scalar> & T,
			const std::size_t componentIndex) const noexcept override;

	std::valarray<scalar> Sv(
			const std::vector<const std::valarray<scalar>*> & concentrations,
			const std::valarray<scalar> & T) const noexcept override;

	std::valarray<scalar> Svk(const std::valarray<scalar> & concentration,
			const std::valarray<scalar> & T,
			const std::size_t componentIndex) const noexcept override;

	/**** For point ****/
	scalar Cv(const std::valarray<scalar> & concentrations) const noexcept
			override;

	scalar pFromUv(const std::valarray<scalar> & concentrations,
			const scalar Uv) const noexcept override;

	scalar UvFromp(const std::valarray<scalar> & concentrations,
			const scalar p) const noexcept override;

	scalar pcFromT(const std::valarray<scalar>&, const scalar T) const noexcept
			override;

	scalar pcFromTk(const scalar, const scalar T,
			const std::size_t) const noexcept override;

	scalar UvcFromT(const std::valarray<scalar> & concentrations,
			const scalar T) const noexcept override;

	scalar UvcFromTk(const scalar, const scalar T,
			const std::size_t componentIndex) const noexcept override;

	scalar TFromUv(const std::valarray<scalar> & concentrations,
			const scalar Uv) const noexcept override;

	scalar dpdrho(const std::valarray<scalar>&, const scalar) const noexcept
			override;

	scalar dpdUv(const std::valarray<scalar> & concentrations,
			const scalar) const noexcept override;

	scalar nonIdeality(const std::valarray<scalar>&,
			const scalar) const noexcept override;

	scalar sqSonicSpeed(const std::valarray<scalar> & concentrations,
			const scalar density, const scalar Uv,
			const scalar pressure) const noexcept override;

	scalar Cp(const std::valarray<scalar> & concentrations,
			const scalar) const noexcept override;

	scalar Cpk(const scalar, const scalar,
			const std::size_t componentIndex) const noexcept override;

	scalar hT(const std::valarray<scalar> & concentrations,
			const scalar T) const noexcept override;

	scalar hkT(const scalar concentration, const scalar T,
			const std::size_t componentIndex) const noexcept override;

	scalar Fv(const std::valarray<scalar> & concentrations,
			const scalar T) const noexcept override;

	scalar Fvk(const scalar concentration, const scalar T,
			const std::size_t componentIndex) const noexcept override;

	scalar Sv(const std::valarray<scalar> & concentrations,
			const scalar T) const noexcept override;

	scalar Svk(const scalar concentration, const scalar T,
			const std::size_t componentIndex) const noexcept override;
};
}  // namespace schemi

#endif /* MIXTUREIDEAL_HPP_ */
