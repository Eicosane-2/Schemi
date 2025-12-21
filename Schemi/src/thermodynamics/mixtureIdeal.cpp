/*
 * mixtureIdeal.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "mixtureIdeal.hpp"

#include "globalConstants.hpp"

schemi::mixtureIdeal::mixtureIdeal() noexcept :
		abstractMixtureThermodynamics(0, 0), M(0), CvArr(0), molecMass(0)
{
}

schemi::mixtureIdeal::mixtureIdeal(const scalar Rin, const scalar hPin,
		const std::valarray<scalar> & Min,
		const std::valarray<scalar> & Cvin) noexcept :
		abstractMixtureThermodynamics(Rin, hPin), M(Min), CvArr(Cvin), molecMass(
				M / NAvogardro)
{
}

schemi::scalar schemi::mixtureIdeal::Rv() const noexcept
{
	return R;
}

const std::valarray<schemi::scalar>& schemi::mixtureIdeal::Mv() const noexcept
{
	return M;
}

const std::valarray<schemi::scalar>& schemi::mixtureIdeal::Cvv() const noexcept
{
	return CvArr;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::Cv(
		const std::vector<const std::valarray<scalar>*> & concentrations) const noexcept
{
	std::valarray<scalar> CvMixtureOutput(0., concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixtureOutput += X[k] * CvArr[k];

	return CvMixtureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::pFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const noexcept
{
	std::valarray<scalar> pressureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., pressureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	for (std::size_t i = 0; i < pressureOutput.size(); ++i)
		pressureOutput[i] = idealFluid::pFromUv(R / CvMixture[i], Uv[i]);

	return pressureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::UvFromp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & p) const noexcept
{
	std::valarray<scalar> UvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., UvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	for (std::size_t i = 0; i < UvOutput.size(); ++i)
		UvOutput[i] = idealFluid::UvFromp(R / CvMixture[i], p[i]);

	return UvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::pcFromT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> pressureConcentrationRatioOutput(
			concentrations[0]->size());

	for (std::size_t i = 0; i < pressureConcentrationRatioOutput.size(); ++i)
		pressureConcentrationRatioOutput[i] = idealFluid::pcFromT(R, T[i]);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::pcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T, const std::size_t) const noexcept
{
	std::valarray<scalar> pressureConcentrationRatioOutput(
			concentration.size());

	for (std::size_t i = 0; i < pressureConcentrationRatioOutput.size(); ++i)
		pressureConcentrationRatioOutput[i] = idealFluid::pcFromT(R, T[i]);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::UvcFromT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> InternalEnergyConcentrationRatioOutput(
			concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0.,
			InternalEnergyConcentrationRatioOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] = idealFluid::UvcFromT(
				CvMixture[i], T[i]);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::UvcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> InternalEnergyConcentrationRatioOutput(
			concentration.size());

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] = idealFluid::UvcFromT(
				CvArr[componentIndex], T[i]);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::TFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const noexcept
{
	std::valarray<scalar> TemperatureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., TemperatureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	for (std::size_t i = 0; i < TemperatureOutput.size(); ++i)
		TemperatureOutput[i] = idealFluid::TFromUv(CvMixture[i], Uv[i],
				(*concentrations[0])[i]);

	return TemperatureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::cFrompT(
		const std::vector<const std::valarray<scalar>*>&,
		const std::valarray<scalar> & p,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> concentratonOutput((p.size()));

	for (std::size_t i = 0; i < concentratonOutput.size(); ++i)
		concentratonOutput[i] = idealFluid::cFrompT(R, p[i], T[i]);

	return concentratonOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::cFrompTk(
		const std::valarray<scalar> & p, const std::valarray<scalar> & T,
		const std::size_t) const noexcept
{
	std::valarray<scalar> concentratonOutput((p.size()));

	for (std::size_t i = 0; i < concentratonOutput.size(); ++i)
		concentratonOutput[i] = idealFluid::cFrompT(R, p[i], T[i]);

	return concentratonOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::dpdrho(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> dpdrhoOutput(concentrations[0]->size());

	for (std::size_t i = 0; i < dpdrhoOutput.size(); ++i)
		dpdrhoOutput[i] = idealFluid::dpdrho();

	return dpdrhoOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::dpdUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> dpdUvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., dpdUvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	for (std::size_t i = 0; i < dpdUvOutput.size(); ++i)
		dpdUvOutput[i] = idealFluid::dpdUv(R / CvMixture[i]);

	return dpdUvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::nonIdeality(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> nonIdealityOutput(concentrations[0]->size());

	for (std::size_t i = 0; i < nonIdealityOutput.size(); ++i)
		nonIdealityOutput[i] = idealFluid::nonIdeality();

	return nonIdealityOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::sqSonicSpeed(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & density, const std::valarray<scalar> & Uv,
		const std::valarray<scalar> & pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::Cp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> CpMixtureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., concentrations[0]->size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	for (std::size_t i = 0; i < CpMixtureOutput.size(); ++i)
		CpMixtureOutput[i] = idealFluid::Cp(CvMixture[i], R);

	return CpMixtureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::Cpk(
		const std::valarray<scalar>&, const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> CpkOutput(T.size());

	for (std::size_t i = 0; i < CpkOutput.size(); ++i)
		CpkOutput[i] = idealFluid::Cp(CvArr[componentIndex], R);

	return CpkOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::hT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::hkT(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::Fv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> FvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = idealFluid::Fmx(hPlanck, molecMass, kB, T[i], rearX[i],
				molecVolumeMixture[i]);

	return FvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::Fvk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> FvOutput(concentration.size());

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = idealFluid::Fv(concentration[i], T[i], R,
				M[componentIndex], hPlanck);

	return FvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::Sv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> SvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = idealFluid::Smx(hPlanck, molecMass, kB, T[i], rearX[i],
				molecVolumeMixture[i]);

	return SvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureIdeal::Svk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> SvOutput(concentration.size());

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = idealFluid::Sv(concentration[i], T[i], R,
				M[componentIndex], hPlanck);

	return SvOutput;
}

schemi::scalar schemi::mixtureIdeal::Cv(
		const std::valarray<scalar> & concentrations) const noexcept
{
	scalar CvMixtureOutput { 0 };

	const auto X = calcMolarFrac(concentrations);

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixtureOutput += X[k] * CvArr[k];

	return CvMixtureOutput;
}

schemi::scalar schemi::mixtureIdeal::pFromUv(
		const std::valarray<scalar> & concentrations,
		const scalar Uv) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	return idealFluid::pFromUv(R / CvMixture, Uv);
}

schemi::scalar schemi::mixtureIdeal::UvFromp(
		const std::valarray<scalar> & concentrations,
		const scalar p) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	return idealFluid::UvFromp(R / CvMixture, p);
}

schemi::scalar schemi::mixtureIdeal::pcFromT(const std::valarray<scalar>&,
		const scalar T) const noexcept
{
	return idealFluid::pcFromT(R, T);
}

schemi::scalar schemi::mixtureIdeal::pcFromTk(const scalar, const scalar T,
		const std::size_t) const noexcept
{
	return idealFluid::pcFromT(R, T);
}

schemi::scalar schemi::mixtureIdeal::UvcFromT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	return idealFluid::UvcFromT(CvMixture, T);
}

schemi::scalar schemi::mixtureIdeal::UvcFromTk(const scalar, const scalar T,
		const std::size_t componentIndex) const noexcept
{
	return idealFluid::UvcFromT(CvArr[componentIndex], T);
}

schemi::scalar schemi::mixtureIdeal::TFromUv(
		const std::valarray<scalar> & concentrations,
		const scalar Uv) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	return idealFluid::TFromUv(CvMixture, Uv, concentrations[0]);
}

schemi::scalar schemi::mixtureIdeal::cFrompT(const std::valarray<scalar>&,
		const scalar p, const scalar T) const noexcept
{
	return idealFluid::cFrompT(R, p, T);
}

schemi::scalar schemi::mixtureIdeal::cFrompTk(const scalar p, const scalar T,
		const std::size_t) const noexcept
{
	return idealFluid::cFrompT(R, p, T);
}

schemi::scalar schemi::mixtureIdeal::dpdrho(const std::valarray<scalar>&,
		const scalar) const noexcept
{
	return idealFluid::dpdrho();
}

schemi::scalar schemi::mixtureIdeal::dpdUv(
		const std::valarray<scalar> & concentrations,
		const scalar) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	return idealFluid::dpdUv(R / CvMixture);
}

schemi::scalar schemi::mixtureIdeal::nonIdeality(const std::valarray<scalar>&,
		const scalar) const noexcept
{
	return idealFluid::nonIdeality();
}

schemi::scalar schemi::mixtureIdeal::sqSonicSpeed(
		const std::valarray<scalar> & concentrations, const scalar density,
		const scalar Uv, const scalar pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

schemi::scalar schemi::mixtureIdeal::Cp(
		const std::valarray<scalar> & concentrations,
		const scalar) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	return idealFluid::Cp(CvMixture, R);
}

schemi::scalar schemi::mixtureIdeal::Cpk(const scalar, const scalar,
		const std::size_t componentIndex) const noexcept
{
	return idealFluid::Cp(CvArr[componentIndex], R);
}

schemi::scalar schemi::mixtureIdeal::hT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

schemi::scalar schemi::mixtureIdeal::hkT(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

schemi::scalar schemi::mixtureIdeal::Fv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	return idealFluid::Fmx(hPlanck, molecMass, kB, T, X, molecVolumeMixture)
			* molecConc;
}

schemi::scalar schemi::mixtureIdeal::Fvk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	return idealFluid::Fv(concentration, T, R, M[componentIndex], hPlanck);
}

schemi::scalar schemi::mixtureIdeal::Sv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	return idealFluid::Smx(hPlanck, molecMass, kB, T, X, molecVolumeMixture)
			* molecConc;
}

schemi::scalar schemi::mixtureIdeal::Svk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	return idealFluid::Sv(concentration, T, R, M[componentIndex], hPlanck);
}
