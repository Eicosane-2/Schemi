/*
 * idealFluid.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "idealFluid.hpp"

#include "intExpPow.hpp"
#include "globalConstants.hpp"

schemi::scalar schemi::idealFluid::pFromUv(const scalar gamma1,
		const scalar Uv) const noexcept
{
	return gamma1 * Uv;
}

schemi::scalar schemi::idealFluid::UvFromp(const scalar gamma1,
		const scalar p) const noexcept
{
	return p / gamma1;
}

schemi::scalar schemi::idealFluid::pcFromT(const scalar R,
		const scalar T) const noexcept
{
	return R * T;
}

schemi::scalar schemi::idealFluid::UvcFromT(const scalar Cv,
		const scalar T) const noexcept
{
	return Cv * T;
}

schemi::scalar schemi::idealFluid::TFromUv(const scalar Cv, const scalar Uv,
		const scalar c) const noexcept
{
	return Uv / (c * Cv);
}

schemi::scalar schemi::idealFluid::dpdrho() const noexcept
{
	return 0;
}

schemi::scalar schemi::idealFluid::dpdUv(const scalar gamma1) const noexcept
{
	return gamma1;
}

schemi::scalar schemi::idealFluid::nonIdeality() const noexcept
{
	return 0;
}

schemi::scalar schemi::idealFluid::Cp(const scalar Cv,
		const scalar R) const noexcept
{
	return Cv + R;
}

schemi::scalar schemi::idealFluid::Fv(const scalar c, const scalar T,
		const scalar R, const scalar M, const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * M * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return -c * R * T * (std::log(nQ / (c * NAvogardro)) + 1);
}

schemi::scalar schemi::idealFluid::Sv(const scalar c, const scalar T,
		const scalar R, const scalar M, const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * M * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return c * R * (log(nQ / (c * NAvogardro)) + 5. / 2.);
}

schemi::scalar schemi::idealFluid::Fmx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	return kB * T
			* ((x * std::log(x + stabilizator)).sum() - 1
					- (x * std::log(Vm / deBroigleLambda3)).sum());
}

schemi::scalar schemi::idealFluid::Smx(const scalar h,
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
