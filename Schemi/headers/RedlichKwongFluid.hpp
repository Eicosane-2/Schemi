/*
 * RedlichKwongFluid.hpp
 *
 *  Created on: 2020/04/18
 *      Author: Maxim Boldyrev
 *
 *      Redlich-Kwong fluid equation of state class.
 */

#ifndef REDLICHKWONGFLUID_HPP_
#define REDLICHKWONGFLUID_HPP_

#include <valarray>

#include "scalar.hpp"

namespace schemi
{
class RedlichKwongFluid
{
	scalar dPdC(const scalar c, const scalar T, const scalar a, const scalar b,
			const scalar R) const noexcept;
	scalar dPdT(const scalar c, const scalar T, const scalar a, const scalar b,
			const scalar R) const noexcept;
	scalar dUvdC(const scalar c, const scalar Cv, const scalar T,
			const scalar a, const scalar b) const noexcept;
	scalar dUvdT(const scalar c, const scalar Cv, const scalar T,
			const scalar a, const scalar b) const noexcept;

public:
	scalar pFromUv(const scalar R, const scalar Cv, const scalar Uv,
			const scalar c, const scalar a, const scalar b) const;

	scalar UvFromp(const scalar R, const scalar Cv, const scalar p,
			const scalar c, const scalar a, const scalar b) const;

	scalar pcFromT(const scalar R, const scalar T, const scalar c,
			const scalar a, const scalar b) const noexcept;

	scalar UvcFromT(const scalar Cv, const scalar T, const scalar c,
			const scalar a, const scalar b) const noexcept;

	scalar TFromUv(const scalar Cv, const scalar Uv, const scalar c,
			const scalar a, const scalar b) const;

	scalar cFrompT(const scalar R, const scalar p, const scalar T,
			const scalar a, const scalar b) const;

	scalar dpdrho(const scalar Cv, const scalar Uv, const scalar c,
			const scalar M, const scalar a, const scalar b,
			const scalar R) const;

	scalar dpdUv(const scalar Cv, const scalar c, const scalar Uv,
			const scalar a, const scalar b, const scalar R) const;

	scalar nonIdeality(const scalar c, const scalar a, const scalar b,
			const scalar T) const noexcept;

	scalar Cp(const scalar Cv, const scalar R, const scalar a, const scalar b,
			const scalar c, const scalar T) const noexcept;

	scalar Fv(const scalar c, const scalar T, const scalar R, const scalar M,
			const scalar a, const scalar b, const scalar h) const noexcept;

	scalar Sv(const scalar c, const scalar T, const scalar R, const scalar M,
			const scalar a, const scalar b, const scalar h) const noexcept;

	scalar Fmx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm, const scalar am, const scalar bm) const noexcept;

	scalar Smx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm, const scalar am, const scalar bm) const noexcept;
};
}  // namespace schemi

#endif /* REDLICHKWONGFLUID_HPP_ */
