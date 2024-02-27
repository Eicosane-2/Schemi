/*
 * KataokaVanDerWaalsFluid.hpp
 *
 *  Created on: 2024/02/18
 *      Author: Maxim Boldyrev
 */

#ifndef KATAOKAVANDERWAALSFLUID_HPP_
#define KATAOKAVANDERWAALSFLUID_HPP_

#include "scalar.hpp"

#include <valarray>

namespace schemi
{
class KataokaVanDerWaalsFluid
{
	constexpr static scalar M { 13 };

	scalar qa(const scalar c, const scalar T, const scalar V0, const scalar eps,
			const scalar R) const noexcept;

	scalar UvA(const scalar c, const scalar T, const scalar V0,
			const scalar eps, const scalar R) const noexcept;

	scalar pA(const scalar c, const scalar T, const scalar V0, const scalar eps,
			const scalar R) const noexcept;

	scalar dPdC(const scalar c, const scalar T, const scalar V0,
			const scalar eps, const scalar R, const scalar b) const noexcept;
	scalar dPdT(const scalar c, const scalar T, const scalar V0,
			const scalar eps, const scalar R, const scalar b) const noexcept;
	scalar dUvdC(const scalar c, const scalar Cv, const scalar T,
			const scalar V0, const scalar eps, const scalar R) const noexcept;
	scalar dUvdT(const scalar c, const scalar Cv, const scalar T,
			const scalar V0, const scalar eps, const scalar R) const noexcept;

public:
	scalar pFromUv(const scalar c, const scalar Cv, const scalar Uv,
			const scalar V0, const scalar eps, const scalar R,
			const scalar b) const;

	scalar UvFromp(const scalar c, const scalar Cv, const scalar p,
			const scalar V0, const scalar eps, const scalar R,
			const scalar b) const;

	scalar pcFromT(const scalar c, const scalar T, const scalar V0,
			const scalar eps, const scalar R, const scalar b) const noexcept;

	scalar UvcFromT(const scalar c, const scalar Cv, const scalar T,
			const scalar V0, const scalar eps, const scalar R) const noexcept;

	scalar TFromUv(const scalar c, const scalar Cv, const scalar Uv,
			const scalar V0, const scalar eps, const scalar R) const;

	scalar cFrompT(const scalar p, const scalar T, const scalar V0,
			const scalar eps, const scalar R, const scalar b) const;

	scalar dpdrho(const scalar MolMass, const scalar c, const scalar Cv,
			const scalar Uv, const scalar V0, const scalar eps, const scalar R,
			const scalar b) const;

	scalar dpdUv(const scalar c, const scalar Cv, const scalar Uv,
			const scalar V0, const scalar eps, const scalar R,
			const scalar b) const;

	scalar nonIdeality(const scalar c, const scalar T, const scalar V0,
			const scalar eps, const scalar R) const noexcept;

	scalar Cp(const scalar c, const scalar T, const scalar V0, const scalar eps,
			const scalar b, const scalar R, const scalar Cv) const noexcept;

	scalar Fv(const scalar c, const scalar T, const scalar R,
			const scalar MolMass, const scalar V0, const scalar eps,
			const scalar b, const scalar h) const noexcept;

	scalar Sv(const scalar c, const scalar T, const scalar R,
			const scalar MolMass, const scalar V0, const scalar eps,
			const scalar b, const scalar h) const noexcept;

	scalar Fmx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm, const scalar V0m, const scalar epsm,
			const scalar bm) const noexcept;

	scalar Smx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm, const scalar V0m, const scalar epsm,
			const scalar bm) const noexcept;
};
}  // namespace schemi

#endif /* KATAOKAVANDERWAALSFLUID_HPP_ */
