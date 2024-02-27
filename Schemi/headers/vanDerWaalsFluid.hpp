/*
 * vanDerWaalsFluid.hpp
 *
 *  Created on: 2019/12/01
 *      Author: Maxim Boldyrev
 *
 *      van der Waals fluid equation of state class.
 */

#ifndef VANDERWAALSFLUID_HPP_
#define VANDERWAALSFLUID_HPP_

#include <valarray>

#include "scalar.hpp"

namespace schemi
{
class vanDerWaalsFluid
{
public:
	scalar pFromUv(const scalar gamma1, const scalar Uv, const scalar c,
			const scalar a, const scalar b) const noexcept;

	scalar UvFromp(const scalar gamma1, const scalar p, const scalar c,
			const scalar a, const scalar b) const noexcept;

	scalar pcFromT(const scalar R, const scalar T, const scalar c,
			const scalar a, const scalar b) const noexcept;

	scalar UvcFromT(const scalar Cv, const scalar T, const scalar c,
			const scalar a) const noexcept;

	scalar TFromUv(const scalar Cv, const scalar Uv, const scalar c,
			const scalar a) const noexcept;

	scalar cFrompT(const scalar R, const scalar p, const scalar T,
			const scalar a, const scalar b) const;

	scalar dpdrho(const scalar gamma1, const scalar Uv, const scalar c,
			const scalar M, const scalar a, const scalar b) const noexcept;

	scalar dpdUv(const scalar gamma1, const scalar c,
			const scalar b) const noexcept;

	scalar nonIdeality(const scalar c, const scalar a) const noexcept;

	scalar Cp(const scalar Cv, const scalar R, const scalar a, const scalar b,
			const scalar c, const scalar T) const noexcept;

	scalar Fv(const scalar c, const scalar T, const scalar R, const scalar M,
			const scalar a, const scalar b, const scalar h) const noexcept;

	scalar Sv(const scalar c, const scalar T, const scalar R, const scalar M,
			const scalar b, const scalar h) const noexcept;

	scalar Fmx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm, const scalar am, const scalar bm) const noexcept;

	scalar Smx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm, const scalar bm) const noexcept;
};
}  // namespace schemi

#endif /* VANDERWAALSFLUID_HPP_ */
