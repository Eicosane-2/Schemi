/*
 * stiffenedFluid.hpp
 *
 *  Created on: 2023/12/07
 *      Author: Maxim Boldyrev
 */

#ifndef STIFFENEDFLUID_HPP_
#define STIFFENEDFLUID_HPP_

#include <valarray>

#include "scalar.hpp"

namespace schemi
{
class stiffenedFluid
{
public:
	scalar pFromUv(const scalar gamma1, const scalar Uv, const scalar p0,
			const scalar gamma) const noexcept;

	scalar UvFromp(const scalar gamma1, const scalar p, const scalar p0,
			const scalar gamma) const noexcept;

	scalar pcFromT(const scalar R, const scalar T, const scalar p0,
			const scalar c) const noexcept;

	scalar UvcFromT(const scalar Cv, const scalar T, const scalar p0,
			const scalar c) const noexcept;

	scalar TFromUv(const scalar Cv, const scalar Uv, const scalar c,
			const scalar p0) const noexcept;

	scalar dpdrho() const noexcept;

	scalar dpdUv(const scalar gamma1) const noexcept;

	scalar nonIdeality(const scalar p0) const noexcept;

	scalar Cp(const scalar Cv, const scalar gamma) const noexcept;

	scalar Fv(const scalar c, const scalar T, const scalar R, const scalar M,
			const scalar p0, const scalar h) const noexcept;

	scalar Sv(const scalar c, const scalar T, const scalar R, const scalar M,
			const scalar h) const noexcept;

	scalar Fmx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm, const scalar p0) const noexcept;

	scalar Smx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm) const noexcept;
};
}  // namespace schemi

#endif /* STIFFENEDFLUID_HPP_ */
