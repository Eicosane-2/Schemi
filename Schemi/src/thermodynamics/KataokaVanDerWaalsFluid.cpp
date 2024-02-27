/*
 * KataokaVanDerWaalsFluid.cpp
 *
 *  Created on: 2024/02/18
 *      Author: Maxim Boldyrev
 */

#include "KataokaVanDerWaalsFluid.hpp"

#include <functional>

#include "intExpPow.hpp"
#include "globalConstants.hpp"
#include "NewtonMethod.hpp"
#include "secantMethod.hpp"

schemi::scalar schemi::KataokaVanDerWaalsFluid::qa(const scalar c,
		const scalar T, const scalar V0, const scalar eps,
		const scalar R) const noexcept
{
	const auto V0c { V0 * c };
	const auto epsRT { 0.5 * eps / (R * T) };

	scalar series { 1.0 };

	for (std::size_t J = 2; J <= M; ++J)
		series += std::pow(V0c, (J - 1.) / 3.) * std::exp((J - 1.) * epsRT);

	return series;
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::UvA(const scalar c,
		const scalar T, const scalar V0, const scalar eps,
		const scalar R) const noexcept
{
	const auto V0c { V0 * c };
	const auto epsRT { 0.5 * eps / (R * T) };

	scalar series { 0 };

	for (std::size_t J = 2; J <= M; ++J)
		series -= (J - 1.) * std::pow(V0c, (J - 1.) / 3.)
				* std::exp((J - 1.) * epsRT);

	return 0.5 * c * eps * series / qa(c, T, V0, eps, R);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::pA(const scalar c,
		const scalar T, const scalar V0, const scalar eps,
		const scalar R) const noexcept
{
	const auto V0c { V0 * c };
	const auto epsRT { 0.5 * eps / (R * T) };

	scalar series { 0 };

	for (std::size_t J = 2; J <= M; ++J)
		series -= (J - 1.) / 3. * std::pow(V0c, (J - 1.) / 3.)
				* std::exp((J - 1.) * epsRT);

	return c * R * T * series / qa(c, T, V0, eps, R);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::dPdC(const scalar c,
		const scalar T, const scalar V0, const scalar eps, const scalar R,
		const scalar b) const noexcept
{
	const auto V0c { V0 * c };
	const auto epsRT { 0.5 * eps / (R * T) };

	scalar numerator13 { 0 };

	for (std::size_t J = 2; M < J; ++J)
		numerator13 -= (J - 1.) / 3. * std::pow(V0c, (J - 1.) / 3.)
				* std::exp((J - 1) * epsRT);

	const auto qa1 = qa(c, T, V0, eps, R);

	return R * T / (1 - c * b) + c * R * T * b / pow<scalar, 2>(1 - c * b)

	+ (R * T * numerator13

	+ V0c * R * T * ([&]()
	{
		scalar series { 0 };

		for (std::size_t J = 2; M < J; ++J)
			series -= pow<scalar, 2>((J - 1.) / 3.) * std::exp((J - 1.) * epsRT)* std::pow(V0c, (J - 4.) / 3.);

		return series;
	}()
	- numerator13 * [&]()
	{
		scalar series
		{	0};

		for (std::size_t J = 2; M < J; ++J)
		series += (J - 1) / 3. * std::pow(V0c, (J - 4.) / 3.) * std::exp((J - 1.) * epsRT);

		return series;
	}() / qa1)
	) / qa1;
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::dPdT(const scalar c,
		const scalar T, const scalar V0, const scalar eps, const scalar R,
		const scalar b) const noexcept
{
	const auto V0c { V0 * c };
	const auto epsRT { 0.5 * eps / (R * T) };

	scalar numerator13 { 0 };

	for (std::size_t J = 2; M < J; ++J)
		numerator13 -= (J - 1.) / 3. * std::pow(V0c, (J - 1.) / 3.)
				* std::exp((J - 1) * epsRT);

	const auto qa1 = qa(c, T, V0, eps, R);

	return c * R / (1 - c * b)

	+ (c * R * numerator13

	+ 0.5 * c * eps / T * ([&]()
	{
		scalar series { 0 };

		for (std::size_t J = 2; M < J; ++J)
			series += pow<scalar, 2>(J - 1.) / 3. * std::pow(V0c, (J - 1.) / 3.)* std::exp((J - 1.) * epsRT);

		return series;
	}()
	- numerator13 * [&]()
	{
		scalar series
		{	0};

		for (std::size_t J = 2; M < J; ++J)
		series -= (J - 1.) * std::pow(V0c, (J - 1.) / 3.) * std::exp((J - 1.) * epsRT);

		return series;
	}() / qa1)
	) / qa1;
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::dUvdC(const scalar c,
		const scalar Cv, const scalar T, const scalar V0, const scalar eps,
		const scalar R) const noexcept
{
	const auto V0c { V0 * c };
	const auto epsRT { 0.5 * eps / (R * T) };

	scalar numerator13 { 0 };

	for (std::size_t J = 2; M < J; ++J)
		numerator13 -= (J - 1.) * std::pow(V0c, (J - 1.) / 3.)
				* std::exp((J - 1) * epsRT);

	const auto qa1 = qa(c, T, V0, eps, R);

	return Cv * T

	+ (0.5 * eps * numerator13

	+ 0.5 * V0c * eps * ([&]()
	{
		scalar series { 0 };

		for (std::size_t J = 2; M < J; ++J)
			series -= pow<scalar, 2>(J - 1.) / 3. * std::pow(V0c, (J - 4.) / 3.)* std::exp((J - 1) * epsRT);

		return series;
	}()
	- numerator13 * [&]()
	{
		scalar series
		{	0};

		for (std::size_t J = 2; M < J; ++J)
		series += (J - 1.) / 3. * std::pow(V0c, (J - 4.) /3. ) * std::exp((J - 1) * epsRT);

		return series;
	}() / qa1)
	) / qa1;
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::dUvdT(const scalar c,
		const scalar Cv, const scalar T, const scalar V0, const scalar eps,
		const scalar R) const noexcept
{
	const auto V0c { V0 * c };
	const auto epsRT { 0.5 * eps / (R * T) };

	const auto qa1 = qa(c, T, V0, eps, R);

	return c * Cv

	+ 0.25 * c * pow<scalar, 2>(eps) / (R * pow<scalar, 2>(T)) * [&]()
	{
		scalar series { 0 };

		for (std::size_t J = 2; M < J; ++J)
			series += pow<scalar, 2>(J - 1.) * std::pow(V0c, (J - 1.) / 3.) * std::exp((J - 1) * epsRT);

		return series;
	}() / qa1

	- 0.25 * c * pow<scalar,2>(eps) / (R * pow<scalar, 2>(T)) *
	pow<scalar, 2>([&]()
			{
				scalar series
				{	0};

				for (std::size_t J = 2; M < J; ++J)
				series -= (J - 1.) * std::pow(V0c, (J - 1.) / 3.) * std::exp((J - 1) * epsRT);

				return series;
			}()) / pow<scalar,2>(qa1);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::pFromUv(const scalar c,
		const scalar Cv, const scalar Uv, const scalar V0, const scalar eps,
		const scalar R, const scalar b) const noexcept
{
	auto U_T = [this, &c, &Cv, &Uv, &V0, &eps, &R](const scalar T) 
	{
		return c * Cv * T + UvA(c, T, V0, eps, R) - Uv;
	};

	auto dU_dT = [this, &c, &Cv, &V0, &eps, &R](const scalar T) 
	{
		return dUvdT(c, Cv, T, V0, eps, R);
	};

	const scalar initT { Uv / ((c + stabilizator) * Cv) };

	const auto T = NewtonMethod(initT, U_T, dU_dT);

	return c * R * T / (1 - c * b) + pA(c, T, V0, eps, R);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::UvFromp(const scalar c,
		const scalar Cv, const scalar p, const scalar V0, const scalar eps,
		const scalar R, const scalar b) const noexcept
{
	auto p_T = [this, &c, &p, &V0, &eps, &b, &R](const scalar T) 
	{
		return c * R * T / (1 - c * b) + pA(c, T, V0, eps, R) - p;
	};

	auto dp_dT = [this, &c, &V0, &eps, &b, &R](const scalar T) 
	{
		return dPdT(c, T, V0, eps, R, b);
	};

	const scalar initT { p / ((c + stabilizator) * R) };

	const auto T = NewtonMethod(initT, p_T, dp_dT);

	return c * Cv * T + UvA(c, T, V0, eps, R);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::pcFromT(const scalar c,
		const scalar T, const scalar V0, const scalar eps, const scalar R,
		const scalar b) const noexcept
{
	return R * T / (1 - c * b) + pA(c, T, V0, eps, R) / (c + stabilizator);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::UvcFromT(const scalar c,
		const scalar Cv, const scalar T, const scalar V0, const scalar eps,
		const scalar R) const noexcept
{
	return Cv * T + UvA(c, T, V0, eps, R) / (c + stabilizator);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::TFromUv(const scalar c,
		const scalar Cv, const scalar Uv, const scalar V0, const scalar eps,
		const scalar R) const noexcept
{
	auto U_T = [this, &c, &Cv, &Uv, &V0, &eps, &R](const scalar T) 
	{
		return c * Cv * T + UvA(c, T, V0, eps, R) - Uv;
	};

	auto dU_dT = [this, &c, &Cv, &V0, &eps, &R](const scalar T) 
	{
		return dUvdT(c, Cv, T, V0, eps, R);
	};

	const scalar initT { Uv / ((c + stabilizator) * Cv) };

	return NewtonMethod(initT, U_T, dU_dT);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::cFrompT(const scalar p,
		const scalar T, const scalar V0, const scalar eps, const scalar R,
		const scalar b) const
{
	auto c_pT = [this, &p, &T, &V0, &eps, &R, &b](const scalar c) 
	{
		return c * R * T / (1 - c * b) + pA(c, T, V0, eps, R) - p;
	};

	const scalar initC { p / (R * T) }, guessC = { 1.1 * initC
			+ 2 * convergenceToleranceGlobal };

	return secantMethod(initC, guessC, c_pT);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::dpdrho(const scalar MolMass,
		const scalar c, const scalar Cv, const scalar Uv, const scalar V0,
		const scalar eps, const scalar R, const scalar b) const noexcept
{
	auto U_T = [this, &c, &Cv, &Uv, &V0, &eps, &R](const scalar T) 
	{
		return c * Cv * T + UvA(c, T, V0, eps, R) - Uv;
	};

	auto dU_dT = [this, &c, &Cv, &V0, &eps, &R](const scalar T) 
	{
		return dUvdT(c, Cv, T, V0, eps, R);
	};

	const scalar initT { Uv / ((c + stabilizator) * Cv) };

	const auto T = NewtonMethod(initT, U_T, dU_dT);

	return (dPdC(c, T, V0, eps, R, b)
			- dPdT(c, T, V0, eps, R, b) * dUvdC(c, Cv, T, V0, eps, R)
					/ dUvdT(c, Cv, T, V0, eps, R)) / MolMass;
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::dpdUv(const scalar c,
		const scalar Cv, const scalar Uv, const scalar V0, const scalar eps,
		const scalar R, const scalar b) const noexcept
{
	auto U_T = [this, &c, &Cv, &Uv, &V0, &eps, &R](const scalar T) 
	{
		return c * Cv * T + UvA(c, T, V0, eps, R) - Uv;
	};

	auto dU_dT = [this, &c, &Cv, &V0, &eps, &R](const scalar T) 
	{
		return dUvdT(c, Cv, T, V0, eps, R);
	};

	const scalar initT { Uv / ((c + stabilizator) * Cv) };

	const auto T = NewtonMethod(initT, U_T, dU_dT);

	return dPdT(c, T, V0, eps, R, b) / dUvdT(c, Cv, T, V0, eps, R);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::nonIdeality(const scalar c,
		const scalar T, const scalar V0, const scalar eps,
		const scalar R) const noexcept
{
	return UvA(c, T, V0, eps, R);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::Cp(const scalar c,
		const scalar T, const scalar V0, const scalar eps, const scalar b,
		const scalar R, const scalar Cv) const noexcept
{
	const auto Vm { 1 / c };
	const auto V0Vm { V0 / Vm };
	const auto epsRT { 0.5 * eps / (R * T) };

	const auto dPdVm { -R * T / pow<scalar, 2>(Vm - b) + (R * T * ([&]()
	{
		scalar series { 0 };

		for (std::size_t J = 2; M < J; ++J)
			series += pow<scalar, 2>((J - 1.) / 3.) * std::pow(V0, (J - 1.) / 3.)
	/ std::pow(Vm, (J - 1.) / 3. + 1.) * std::exp((J - 1.) * epsRT);

		return series;
	}()))/(Vm * (qa(c, T, V0, eps, R)))
	-(R * T * ([&]()
	{
		scalar series
		{	0};

		for (std::size_t J = 2; M < J; ++J)
		series -= (J - 1.) / 3. * std::pow(V0Vm, (J - 1.) / 3.) * std::exp((J - 1.) * epsRT);

		return series;
	}()))/(pow<scalar, 2>(Vm)*(qa(c, T, V0, eps, R)))
	-(R * T * ([&]()
	{
		scalar series
		{	0};

		for (std::size_t J = 2; M < J; ++J)
		series -= (J - 1.) / 3. * std::pow(V0Vm, (J - 1.) / 3.) * std::exp((J - 1.) * epsRT);

		return series;
	}()) * ([&]()
	{
		scalar series
		{	0};

		for (std::size_t J = 2; M < J; ++J)
		series -= (J - 1.) / 3. * std::pow(V0, (J - 1.) / 3.) / std::pow(Vm, (J - 1.) / 3. + 1.) * std::exp((J - 1.) * epsRT);

		return series;
	}()))/(Vm * (pow<scalar, 2>(qa(c, T, V0, eps, R))))
};

	return Cv - T * pow<scalar, 2>(dPdT(c, T, V0, eps, R, b)) / dPdVm;
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::Fv(const scalar c,
		const scalar T, const scalar R, const scalar MolMass, const scalar V0,
		const scalar eps, const scalar b, const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * MolMass * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return -c * R * T
			* (std::log(nQ * (1. / (c + stabilizator) - b) / NAvogardro) + 1)
			- c * R * T * std::log(qa(c, T, V0, eps, R));
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::Sv(const scalar c,
		const scalar T, const scalar R, const scalar MolMass, const scalar V0,
		const scalar eps, const scalar b, const scalar h) const noexcept
{
	const auto sqrt_nQ { std::sqrt(
			(2 * Pi_number * MolMass * R * T) / pow<scalar, 2>(NAvogardro * h)) };
	const auto nQ { pow<scalar, 3>(sqrt_nQ) };

	return c * R
			* (log(nQ * (1. / (c + stabilizator) - b) / NAvogardro) + 5. / 2.)
			+ c * R * std::log(qa(c, T, V0, eps, R)) + c * eps / T * [&]()
			{
				const auto V0c { V0 * c };
				const auto epsRT { 0.5 * eps / (R * T) };

				scalar series { 0.0 };

				for (std::size_t J = 2; J <= M; ++J)
					series -= 0.5 * (J - 1) * std::pow(V0c, (J - 1) / 3.) * std::exp(0.5 * (J - 1.) * epsRT);

				return series;
			}() / (qa(c, T, V0, eps, R));
		}

schemi::scalar schemi::KataokaVanDerWaalsFluid::Fmx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm, const scalar V0m,
		const scalar epsm, const scalar bm) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	const auto mixF { (x * std::log(x + stabilizator)).sum() };
	const auto transF { (x * std::log((Vm - bm) / deBroigleLambda3)).sum() };
	const auto potF { std::log(qa(1 / Vm, T, V0m, epsm, kB)) };

	return kB * T * ((mixF - 1 - transF) - potF);
}

schemi::scalar schemi::KataokaVanDerWaalsFluid::Smx(const scalar h,
		const std::valarray<scalar> & massx, const scalar kB, const scalar T,
		const std::valarray<scalar> & x, const scalar Vm, const scalar V0m,
		const scalar epsm, const scalar bm) const noexcept
{
	const auto deBroigleLambda = std::sqrt(
			pow<scalar, 2>(h) / (2 * Pi_number * massx * kB * T));
	const auto deBroigleLambda3 = deBroigleLambda * deBroigleLambda
			* deBroigleLambda;

	const auto mixS { -(x * std::log(x + stabilizator)).sum() };
	const auto transS { (x * std::log((Vm - bm) / deBroigleLambda3)).sum() };
	const auto potS { std::log(qa(1 / Vm, T, V0m, epsm, kB))
			+ epsm / (kB * T) * [&]()
			{
				const auto V0c { V0m / Vm };
				const auto epsRT { 0.5 * epsm / (kB * T) };

				scalar series { 0.0 };

				for (std::size_t J = 2; J <= M; ++J)
					series -= 0.5 * (J - 1) * std::pow(V0c, (J - 1) / 3.) * std::exp(0.5 * (J - 1.) * epsRT);

				return series;
			}() / (qa(1 / Vm, T, V0m, epsm, kB))};

	return kB * ((mixS + 5. / 2. + transS) - potS);
}
