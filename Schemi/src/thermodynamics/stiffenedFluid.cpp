/*
 * stiffenedFluid.cpp
 *
 *  Created on: 2023/12/07
 *      Author: Maxim Boldyrev
 */

#include "stiffenedFluid.hpp"

#include "globalConstants.hpp"
#include "intExpPow.hpp"

schemi::scalar schemi::stiffenedFluid::pFromUv(const scalar gamma1,
		const scalar Uv, const scalar p0, const scalar gamma) const noexcept
{
	return gamma1 * Uv - p0 * gamma;
}

schemi::scalar schemi::stiffenedFluid::UvFromp(const scalar gamma1,
		const scalar p, const scalar p0, const scalar gamma) const noexcept
{
	return (p + p0 * gamma) / gamma1;
}

schemi::scalar schemi::stiffenedFluid::pcFromT(const scalar R, const scalar T,
		const scalar p0, const scalar c) const noexcept
{
	return R * T - p0 / c;
}

schemi::scalar schemi::stiffenedFluid::UvcFromT(const scalar Cv, const scalar T,
		const scalar p0, const scalar c) const noexcept
{
	return Cv * T + p0 / c;
}

schemi::scalar schemi::stiffenedFluid::TFromUv(const scalar Cv, const scalar Uv,
		const scalar c, const scalar p0) const noexcept
{
	return (Uv - p0) / (c * Cv);
}

schemi::scalar schemi::stiffenedFluid::dpdrho() const noexcept
{
	return 0;
}

schemi::scalar schemi::stiffenedFluid::dpdUv(const scalar gamma1) const noexcept
{
	return gamma1;
}

schemi::scalar schemi::stiffenedFluid::nonIdeality(
		const scalar p0) const noexcept
{
	return p0;
}

schemi::scalar schemi::stiffenedFluid::Cp(const scalar Cv,
		const scalar gamma) const noexcept
{
	return Cv * gamma;
}

schemi::scalar schemi::stiffenedFluid::Fv(const scalar c, const scalar T,
		const scalar R, const scalar M, const scalar p0,
		const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * M * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return -c * R * T * (std::log(nQ / (c * NAvogardro)) + 1) + p0;
}

schemi::scalar schemi::stiffenedFluid::Sv(const scalar c, const scalar T,
		const scalar R, const scalar M, const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * M * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return c * R * (log(nQ / (c * NAvogardro)) + 5. / 2.);
}

schemi::scalar schemi::stiffenedFluid::Fmx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm,
		const scalar p0) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	return kB * T
			* ((x * std::log(x + stabilizator)).sum() - 1
					- (x * std::log(Vm / deBroigleLambda3)).sum()) + p0 * Vm;
}

schemi::scalar schemi::stiffenedFluid::Smx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	return kB
			* ((-x * std::log(x + stabilizator)).sum()
					+ (x * std::log(Vm / deBroigleLambda3)).sum() + 5. / 2.);
}
