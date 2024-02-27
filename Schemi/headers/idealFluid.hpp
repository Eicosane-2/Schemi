/*
 * idealFluid.hpp
 *
 *  Created on: 2019/12/01
 *      Author: Maxim Boldyrev
 *
 *      Ideal fluid equation of state (Mendeleev-Clapeyron ideal gas) class.
 */

#ifndef IDEALFLUID_HPP_
#define IDEALFLUID_HPP_

#include <valarray>

#include "scalar.hpp"

namespace schemi
{
class idealFluid
{
public:
	scalar pFromUv(const scalar gamma1, const scalar Uv) const noexcept;

	scalar UvFromp(const scalar gamma1, const scalar p) const noexcept;

	scalar pcFromT(const scalar R, const scalar T) const noexcept;

	scalar UvcFromT(const scalar Cv, const scalar T) const noexcept;

	scalar TFromUv(const scalar Cv, const scalar Uv,
			const scalar c) const noexcept;

	scalar cFrompT(const scalar R, const scalar p,
			const scalar T) const noexcept;

	scalar dpdrho() const noexcept;

	scalar dpdUv(const scalar gamma1) const noexcept;

	scalar nonIdeality() const noexcept;

	scalar Cp(const scalar Cv, const scalar R) const noexcept;

	scalar Fv(const scalar c, const scalar T, const scalar R, const scalar M,
			const scalar h) const noexcept;

	scalar Sv(const scalar c, const scalar T, const scalar R, const scalar M,
			const scalar h) const noexcept;

	scalar Fmx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm) const noexcept;

	scalar Smx(const scalar h, const std::valarray<scalar> & massx,
			const scalar kB, const scalar T, const std::valarray<scalar> & x,
			const scalar Vm) const noexcept;
};
}  // namespace schemi

#endif /* IDEALFLUID_HPP_ */
