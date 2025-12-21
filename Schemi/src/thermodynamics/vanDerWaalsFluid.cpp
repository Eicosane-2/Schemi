/*
 * vanDerWaalsFluid.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "vanDerWaalsFluid.hpp"

#include "cubicEquationSolver.hpp"
#include "globalConstants.hpp"
#include "intExpPow.hpp"

schemi::scalar schemi::vanDerWaalsFluid::pFromUv(const scalar gamma1,
		const scalar Uv, const scalar c, const scalar a,
		const scalar b) const noexcept
{
	const auto acc = a * pow<scalar, 2>(c);

	return gamma1 * (Uv + acc) / (1 - c * b) - acc;
}

schemi::scalar schemi::vanDerWaalsFluid::UvFromp(const scalar gamma1,
		const scalar p, const scalar c, const scalar a,
		const scalar b) const noexcept
{
	const auto acc = a * pow<scalar, 2>(c);

	return (p + acc) * (1 - c * b) / gamma1 - acc;
}

schemi::scalar schemi::vanDerWaalsFluid::pcFromT(const scalar R, const scalar T,
		const scalar c, const scalar a, const scalar b) const noexcept
{
	return R * T / (1 - c * b) - a * c;
}

schemi::scalar schemi::vanDerWaalsFluid::UvcFromT(const scalar Cv,
		const scalar T, const scalar c, const scalar a) const noexcept
{
	return Cv * T - a * c;
}

schemi::scalar schemi::vanDerWaalsFluid::TFromUv(const scalar Cv,
		const scalar Uv, const scalar c, const scalar a) const noexcept
{
	return (Uv + a * pow<scalar, 2>(c)) / ((c + stabilizator) * Cv);
}

schemi::scalar schemi::vanDerWaalsFluid::cFrompT(const scalar R, const scalar p,
		const scalar T, const scalar a, const scalar b) const
{
	return returnSinglePosValue(
			cubicEquationSolver(a * b, -a, R * T + b * p, -p));
}

schemi::scalar schemi::vanDerWaalsFluid::dpdrho(const scalar gamma1,
		const scalar Uv, const scalar c, const scalar M, const scalar a,
		const scalar b) const noexcept
{
	const auto m_cb_p_1 { 1 - c * b };

	const auto ac { a * c };

	const auto dpdrho { (gamma1
			* (2 * ac / m_cb_p_1
					- b * (Uv + ac * c) / (pow<scalar, 2>(m_cb_p_1))) - 2 * ac)
			/ M };
	return dpdrho;
}

schemi::scalar schemi::vanDerWaalsFluid::dpdUv(const scalar gamma1,
		const scalar c, const scalar b) const noexcept
{
	return gamma1 / (1 - c * b);
}

schemi::scalar schemi::vanDerWaalsFluid::nonIdeality(const scalar c,
		const scalar a) const noexcept
{
	return -a * pow<scalar, 2>(c);
}

schemi::scalar schemi::vanDerWaalsFluid::Cp(const scalar Cv, const scalar R,
		const scalar a, const scalar b, const scalar c,
		const scalar T) const noexcept
{
	const auto m_cb_p_1 { 1 - c * b };

	const auto R_modif { R
			/ (1 - 2 * a * c * pow<scalar, 2>(m_cb_p_1) / (R * T)) };

	return Cv + R_modif;
}

schemi::scalar schemi::vanDerWaalsFluid::Fv(const scalar c, const scalar T,
		const scalar R, const scalar M, const scalar a, const scalar b,
		const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * M * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return -c * R * T
			* (std::log(nQ * (1. / (c + stabilizator) - b) / NAvogardro) + 1)
			- a * c * c;
}

schemi::scalar schemi::vanDerWaalsFluid::Sv(const scalar c, const scalar T,
		const scalar R, const scalar M, const scalar b,
		const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * M * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return c * R
			* (std::log(nQ * (1. / (c + stabilizator) - b) / NAvogardro)
					+ 5. / 2.);
}

schemi::scalar schemi::vanDerWaalsFluid::Fmx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm, const scalar am,
		const scalar bm) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	return kB * T
			* ((x * std::log(x + stabilizator)).sum() - 1
					- (x * std::log((Vm - bm) / deBroigleLambda3)).sum())
			- am / Vm;
}

schemi::scalar schemi::vanDerWaalsFluid::Smx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm,
		const scalar bm) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	return kB
			* ((-x * std::log(x + stabilizator)).sum() + 5. / 2.
					+ (x * std::log((Vm - bm) / deBroigleLambda3)).sum());
}
