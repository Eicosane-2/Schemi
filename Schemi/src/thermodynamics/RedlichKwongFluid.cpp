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
	return sqrtTemp * sqrtTemp;
}

schemi::scalar schemi::RedlichKwongFluid::dpdrho(const scalar Cv,
		const scalar gamma1, const scalar Uv, const scalar c, const scalar M,
		const scalar a, const scalar b) const
{
	const auto sqrtT { returnSinglePosValue(
			cubicEquationSolver(c * Cv, 0, -Uv,
					-1.5 * a * c / b * std::log1p(c * b))) };

	const auto logar { std::log1p(c * b) };

	const auto m_cb_p_1 { 1 - c * b };
	const auto cb_p_1 { 1 + c * b };
	const auto ac { a * c };

	return (gamma1 * b / (m_cb_p_1 * m_cb_p_1)
			* (Uv + 1.5 * ac * logar / (sqrtT * b))
			+ gamma1 / m_cb_p_1
					* (1.5 * a * logar / (b * sqrtT)
							+ 0.75 * a * Uv * logar
									/ (c * b * Cv * pow<scalar, 3>(sqrtT))
							+ 1.5 * ac / (sqrtT * cb_p_1))
			- 2 * ac / (sqrtT * cb_p_1)
			- 0.5 * a * Uv / (Cv * pow<scalar, 3>(sqrtT) * cb_p_1)
			+ ac * c * b / (sqrtT * pow<scalar, 2>(cb_p_1))) / M;
}

schemi::scalar schemi::RedlichKwongFluid::dpdUv(const scalar Cv,
		const scalar gamma1, const scalar c, const scalar Uv, const scalar a,
		const scalar b) const
{
	const auto sqrtT { returnSinglePosValue(
			cubicEquationSolver(c * Cv, 0, -Uv,
					-1.5 * a * c / b * std::log1p(c * b))) };

	const auto logar { std::log1p(c * b) };

	return gamma1 / (1 - c * b)
			* (1 - 0.75 * a * logar / (Cv * b * pow<scalar, 3>(sqrtT)))
			+ 0.5 * a * c / (Cv * pow<scalar, 3>(sqrtT) * (1 + c * b));
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

	return -c * R * T * (std::log(nQ * (1. / c - b) / NAvogardro) + 1)
			- a * c / (b * std::sqrt(T)) * std::log1p(c * b);
}

schemi::scalar schemi::RedlichKwongFluid::Sv(const scalar c, const scalar T,
		const scalar R, const scalar M, const scalar a, const scalar b,
		const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * M * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return c * R * (log(nQ * (1. / c - b) / NAvogardro) + 5. / 2.)
			- 0.5 * a * c / (b * T * std::sqrt(T)) * std::log1p(c * b);
}

schemi::scalar schemi::RedlichKwongFluid::Fmx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm,
		const std::valarray<std::valarray<scalar>> & a,
		const std::valarray<std::valarray<scalar>> & b) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	scalar bm { 0 };
	for (std::size_t k = 0; k < x.size(); ++k)
		for (std::size_t l = 0; l < x.size(); ++l)
			bm += x[k] * x[l] * b[k][l];

	const auto mixF { (x * std::log(x + stabilizator)).sum() };
	const auto transF { (x * std::log((Vm - bm) / deBroigleLambda3)).sum() };

	scalar potF { 0 };
	const auto sqrtT { std::sqrt(T) };
	for (std::size_t k = 0; k < x.size(); ++k)
		for (std::size_t l = 0; l < x.size(); ++l)
			potF += x[k] * a[k][l] / b[k][l] * std::log1p(b[k][l] * x[l] / Vm);
	potF /= sqrtT;

	return kB * T * (mixF - 1 - transF) - potF;
}

schemi::scalar schemi::RedlichKwongFluid::Smx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm,
		const std::valarray<std::valarray<scalar>> & a,
		const std::valarray<std::valarray<scalar>> & b) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	scalar bm { 0 };
	for (std::size_t k = 0; k < x.size(); ++k)
		for (std::size_t l = 0; l < x.size(); ++l)
			bm += x[k] * x[l] * b[k][l];

	const auto mixS { -(x * std::log(x + stabilizator)).sum() };
	const auto transS { (x * std::log((Vm - bm) / deBroigleLambda3)).sum() };

	scalar potS { 0 };
	const auto TsqrtT { T * std::sqrt(T) };
	for (std::size_t k = 0; k < x.size(); ++k)
		for (std::size_t l = 0; l < x.size(); ++l)
			potS += x[k] * a[k][l] / b[k][l] * std::log1p(b[k][l] * x[l] / Vm);
	potS *= 0.5 / TsqrtT;

	return kB * (mixS + 5. / 2. + transS) - potS;
}
