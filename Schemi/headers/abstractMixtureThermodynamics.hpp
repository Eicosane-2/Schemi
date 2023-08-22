/*
 * abstractMixtureThermodynamics.hpp
 *
 *  Created on: 2020/02/17
 *      Author: Maxim Boldyrev
 *
 *      Interface class for mixture thermodynamics.
 */

#ifndef ABSTRACTMIXTURETHERMODYNAMICS_HPP_
#define ABSTRACTMIXTURETHERMODYNAMICS_HPP_

#include <valarray>
#include <vector>

#include "fractionCalculation.hpp"
#include "scalar.hpp"

namespace schemi
{
class abstractMixtureThermodynamics: public fractionCalculation
{
protected:
	const scalar R, kB, hPlanck;
public:
	abstractMixtureThermodynamics(const scalar Rin, const scalar hPin) noexcept;

	virtual ~abstractMixtureThermodynamics() noexcept =0;

	/* energy/mole/temperature */
	virtual scalar Rv() const noexcept =0;

	/* mass/mole */
	virtual const std::valarray<scalar>& Mv() const noexcept =0;

	/* energy/(mole*temperature) */
	virtual const std::valarray<scalar>& Cvv() const noexcept =0;

	/**** For field ****/
	/* energy/(mole*temperature) */
	virtual std::valarray<scalar> Cv(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/) const noexcept =0;

	/* pressure */
	virtual std::valarray<scalar> pFromUv(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*Uv*/) const =0;

	/* energy/volume */
	virtual std::valarray<scalar> UvFromp(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*p*/) const =0;

	/* pressure/concentration */
	virtual std::valarray<scalar> pcFromT(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*T*/) const noexcept =0;

	/* pressure/concentration */
	virtual std::valarray<scalar> pcFromTk(
			const std::valarray<scalar>& /*concentration*/,
			const std::valarray<scalar>& /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* energy/mole */
	virtual std::valarray<scalar> UvcFromT(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*T*/) const noexcept =0;

	/* energy/mole */
	virtual std::valarray<scalar> UvcFromTk(
			const std::valarray<scalar>& /*concentration*/,
			const std::valarray<scalar>& /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* temperature */
	virtual std::valarray<scalar> TFromUv(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*Uv*/) const =0;

	/* pressure/density */
	virtual std::valarray<scalar> dpdrho(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*Uv*/) const =0;

	/* dimensionless */
	virtual std::valarray<scalar> dpdUv(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*Uv*/) const =0;

	/* energy/volume */
	virtual std::valarray<scalar> nonIdeality(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*T*/) const noexcept =0;

	/* length^2/time^2 */
	virtual std::valarray<scalar> sqSonicSpeed(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*density*/,
			const std::valarray<scalar>& /*Uv*/,
			const std::valarray<scalar>& /*pressure*/) const noexcept =0;

	/* energy/(mole*temperature) */
	virtual std::valarray<scalar> Cp(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*T*/) const noexcept =0;

	/* energy/(mole*temperature) */
	virtual std::valarray<scalar> Cpk(
			const std::valarray<scalar>& /*concentrations*/,
			const std::valarray<scalar>& /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* energy/(mole*temperature) */
	virtual std::valarray<scalar> hT(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*T*/) const noexcept =0;

	/* energy/(mole*temperature) */
	virtual std::valarray<scalar> hkT(
			const std::valarray<scalar>& /*concentration*/,
			const std::valarray<scalar>& /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* energy/volume */
	virtual std::valarray<scalar> Fv(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*T*/) const noexcept =0;

	/* energy/volume */
	virtual std::valarray<scalar> Fvk(
			const std::valarray<scalar>& /*concentration*/,
			const std::valarray<scalar>& /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* energy/volume/temperature */
	virtual std::valarray<scalar> Sv(
			const std::vector<const std::valarray<scalar>*>& /*concentrations*/,
			const std::valarray<scalar>& /*T*/) const noexcept =0;

	/* energy/volume/temperature */
	virtual std::valarray<scalar> Svk(
			const std::valarray<scalar>& /*concentration*/,
			const std::valarray<scalar>& /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/**** For point ****/
	/* energy/(mole*temperature) */
	virtual scalar Cv(
			const std::valarray<scalar>& /*concentrations*/) const noexcept =0;

	/* pressure */
	virtual scalar pFromUv(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*Uv*/) const =0;

	/* energy/volume */
	virtual scalar UvFromp(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*p*/) const =0;

	/* pressure/concentration */
	virtual scalar pcFromT(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*T*/) const noexcept =0;

	/* pressure/concentration */
	virtual scalar pcFromTk(const scalar /*concentration*/, const scalar /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* energy/mole */
	virtual scalar UvcFromT(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*T*/) const noexcept =0;

	/* energy/mole */
	virtual scalar UvcFromTk(const scalar /*concentration*/, const scalar /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* temperature */
	virtual scalar TFromUv(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*Uv*/) const =0;

	/* pressure/density */
	virtual scalar dpdrho(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*Uv*/) const =0;

	/* dimensionless */
	virtual scalar dpdUv(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*Uv*/) const =0;

	/* energy/volume */
	virtual scalar nonIdeality(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*T*/) const noexcept =0;

	/* length^2/time^2 */
	virtual scalar sqSonicSpeed(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*density*/, const scalar /*Uv*/,
			const scalar /*pressure*/) const noexcept =0;

	/* energy/(mole*temperature) */
	virtual scalar Cp(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*T*/) const noexcept =0;

	/* energy/(mole*temperature) */
	virtual scalar Cpk(const scalar /*concentration*/, const scalar /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* energy/(mole*temperature) */
	virtual scalar hT(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*T*/) const noexcept =0;

	/* energy/(mole*temperature) */
	virtual scalar hkT(const scalar /*concentration*/, const scalar /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* energy/volume */
	virtual scalar Fv(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*T*/) const noexcept =0;

	/* energy/volume */
	virtual scalar Fvk(const scalar /*concentration*/, const scalar /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;

	/* energy/volume/temperature */
	virtual scalar Sv(const std::valarray<scalar>& /*concentrations*/,
			const scalar /*T*/) const noexcept =0;

	/* energy/volume/temperature */
	virtual scalar Svk(const scalar /*concentration*/, const scalar /*T*/,
			const std::size_t /*componentIndex*/) const noexcept =0;
};
}  // namespace schemi

#endif /* ABSTRACTMIXTURETHERMODYNAMICS_HPP_ */
