/*
 * RedlichKwongFluid.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "RedlichKwongFluid.hpp"

#include "cubicEquationSolver.hpp"
#include "globalConstants.hpp"
#include "intExpPow.hpp"

schemi::scalar schemi::RedlichKwongFluid::dPdC(const scalar c, const scalar T,
		const scalar a, const scalar b, const scalar R) const noexcept
{
	return R * T / (1 - c * b) + c * R * T * b / pow<scalar, 2>(1 - c * b)
			- 2 * a * c / (std::sqrt(T) * (1 + b * c))
			+ a * pow<scalar, 2>(c) * b
					/ (std::sqrt(T) * pow<scalar, 2>(1 + b * c));
}

schemi::scalar schemi::RedlichKwongFluid::dPdT(const scalar c, const scalar T,
		const scalar a, const scalar b, const scalar R) const noexcept
{
	return c * R / (1 - c * b)
			+ 0.5 * a * pow<scalar, 2>(c)
					/ (pow<scalar, 3>(std::sqrt(T)) * (1 + c * b));
}

schemi::scalar schemi::RedlichKwongFluid::dUvdC(const scalar c, const scalar Cv,
		const scalar T, const scalar a, const scalar b) const noexcept
{
	return Cv * T - 1.5 * a * std::log1p(c * b) / (b * sqrt(T))
			- 1.5 * a * c / (sqrt(T) * (1 + c * b));
}

schemi::scalar schemi::RedlichKwongFluid::dUvdT(const scalar c, const scalar Cv,
		const scalar T, const scalar a, const scalar b) const noexcept
{
	return c * Cv
			+ 0.75 * a * c * std::log1p(c * b)
					/ (b * pow<scalar, 3>(std::sqrt(T)));
}

schemi::scalar schemi::RedlichKwongFluid::pFromUv(const scalar R,
		const scalar Cv, const scalar Uv, const scalar c, const scalar a,
		const scalar b) const
{
	const auto sqrtTemp { returnSinglePosValue(
			cubicEquationSolver(c * Cv, 0, -Uv,
					-1.5 * a * c / b * std::log1p(c * b))) };

	return c * R * pow<scalar, 2>(sqrtTemp) / (1 - c * b)
			- a * pow<scalar, 2>(c) / (sqrtTemp * (1 + c * b));
}

schemi::scalar schemi::RedlichKwongFluid::UvFromp(const scalar R,
		const scalar Cv, const scalar p, const scalar c, const scalar a,
		const scalar b) const
{
	const auto sqrtTemp { returnSinglePosValue(
			cubicEquationSolver(c * R, 0, -p * (1 - c * b),
					-a * pow<scalar, 2>(c) * (1 - c * b) / (1 + c * b))) };

	return c * Cv * pow<scalar, 2>(sqrtTemp)
			- 1.5 * a * c / (b * sqrtTemp) * std::log1p(c * b);
}

schemi::scalar schemi::RedlichKwongFluid::pcFromT(const scalar R,
		const scalar T, const scalar c, const scalar a,
		const scalar b) const noexcept
{
	return R * T / (1 - c * b) - a * c / (std::sqrt(T) * (1 + c * b));
}

schemi::scalar schemi::RedlichKwongFluid::UvcFromT(const scalar Cv,
		const scalar T, const scalar c, const scalar a,
		const scalar b) const noexcept
{
	return Cv * T - 1.5 * a / (b * std::sqrt(T)) * std::log1p(c * b);
}

schemi::scalar schemi::RedlichKwongFluid::TFromUv(const scalar Cv,
		const scalar Uv, const scalar c, const scalar a, const scalar b) const
{
	const auto sqrtTemp { returnSinglePosValue(
			cubicEquationSolver(c * Cv, 0, -Uv,
					-1.5 * a * c / b * std::log1p(c * b))) };
	return pow<scalar, 2>(sqrtTemp);
}

schemi::scalar schemi::RedlichKwongFluid::cFrompT(const scalar R,
		const scalar p, const scalar T, const scalar a, const scalar b) const
{
	return returnSinglePosValue(
			cubicEquationSolver(a * b / std::sqrt(T),
					p * pow<scalar, 2>(b) + R * T * b - a / std::sqrt(T), R * T,
					-p));
}

schemi::scalar schemi::RedlichKwongFluid::dpdrho(const scalar Cv,
		const scalar Uv, const scalar c, const scalar M, const scalar a,
		const scalar b, const scalar R) const
{
	const auto sqrtT { returnSinglePosValue(
			cubicEquationSolver(c * Cv, 0, -Uv,
					-1.5 * a * c / b * std::log1p(c * b))) };

	const auto T = pow<scalar, 2>(sqrtT);

	return (dPdC(c, T, a, b, R)
			- dPdT(c, T, a, b, R) * dUvdC(c, Cv, T, a, b)
					/ dUvdT(c, Cv, T, a, b)) / M;
}

schemi::scalar schemi::RedlichKwongFluid::dpdUv(const scalar Cv, const scalar c,
		const scalar Uv, const scalar a, const scalar b, const scalar R) const
{
	const auto sqrtT { returnSinglePosValue(
			cubicEquationSolver(c * Cv, 0, -Uv,
					-1.5 * a * c / b * std::log1p(c * b))) };

	const auto T = pow<scalar, 2>(sqrtT);

	return dPdT(c, T, a, b, R) / dUvdT(c, Cv, T, a, b);
}

schemi::scalar schemi::RedlichKwongFluid::nonIdeality(const scalar c,
		const scalar a, const scalar b, const scalar T) const noexcept
{
	return -1.5 * a * c / (b * std::sqrt(T)) * std::log1p(c * b);
}

schemi::scalar schemi::RedlichKwongFluid::Cp(const scalar Cv, const scalar R,
		const scalar a, const scalar b, const scalar c,
		const scalar T) const noexcept
{
	const auto sqrtT = std::sqrt(T);

	const auto m_cb_p_1 { 1 - c * b };
	const auto cb_p_1 { 1 + c * b };

	const auto numerator { R / m_cb_p_1 + 0.5 * a * c / (sqrtT * T * cb_p_1) };

	const auto R_modif { T * numerator * numerator
			/ (R * T / (m_cb_p_1 * m_cb_p_1) - a * c / (cb_p_1 * sqrtT)
					- a * c / (cb_p_1 * cb_p_1 * sqrtT)) };
	return Cv + R_modif;
}

schemi::scalar schemi::RedlichKwongFluid::Fv(const scalar c, const scalar T,
		const scalar R, const scalar M, const scalar a, const scalar b,
		const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * M * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return -c * R * T
			* (std::log(nQ * (1. / (c + stabilizator) - b) / NAvogardro) + 1)
			- a * c / (b * std::sqrt(T)) * std::log1p(c * b);
}

schemi::scalar schemi::RedlichKwongFluid::Sv(const scalar c, const scalar T,
		const scalar R, const scalar M, const scalar a, const scalar b,
		const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * M * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return c * R
			* (log(nQ * (1. / (c + stabilizator) - b) / NAvogardro) + 5. / 2.)
			- 0.5 * a * c / (b * T * std::sqrt(T)) * std::log1p(c * b);
}

schemi::scalar schemi::RedlichKwongFluid::Fmx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm, const scalar am,
		const scalar bm) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	const auto mixF { (x * std::log(x + stabilizator)).sum() };
	const auto transF { (x * std::log((Vm - bm) / deBroigleLambda3)).sum() };
	const auto potF { am / (bm * std::sqrt(T)) * std::log1p(bm / Vm) };

	return kB * T * (mixF - 1 - transF) - potF;
}

schemi::scalar schemi::RedlichKwongFluid::Smx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm, const scalar am,
		const scalar bm) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	const auto mixS { -(x * std::log(x + stabilizator)).sum() };
	const auto transS { (x * std::log((Vm - bm) / deBroigleLambda3)).sum() };
	const auto potS { 0.5 * am / (bm * T * std::sqrt(T)) * std::log1p(bm / Vm) };

	return kB * (mixS + 5. / 2. + transS) - potS;
}
