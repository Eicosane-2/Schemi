/*
 * mixtureStiffened.cpp
 *
 *  Created on: 2023/12/11
 *      Author: Maxim Boldyrev
 */

#include "mixtureStiffened.hpp"

#include "globalConstants.hpp"

schemi::mixtureStiffened::mixtureStiffened() noexcept :
		abstractMixtureThermodynamics(0, 0), M(0), CvArr(0), molecMass(0), p0Matrix(
				0), gammaMatrix(0)
{
}

schemi::mixtureStiffened::mixtureStiffened(const scalar Rin, const scalar hPin,
		const std::valarray<scalar> & Min, const std::valarray<scalar> & Cvin,
		const std::valarray<scalar> & p0in,
		const std::valarray<scalar> & gammain) noexcept :
		abstractMixtureThermodynamics(Rin, hPin), M(Min), CvArr(Cvin), molecMass(
				M / NAvogardro), p0Matrix(std::valarray<scalar>(Min.size()),
				Min.size()), gammaMatrix(std::valarray<scalar>(Min.size()),
				Min.size())
{
	const std::size_t numberOfComponents { Min.size() };

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		p0Matrix[k][k] = p0in[k];
		gammaMatrix[k][k] = gammain[k];
	}

	for (std::size_t k = 0; k < numberOfComponents; ++k)
		for (std::size_t l = 0; l < numberOfComponents; ++l)
			if (k != l)
			{
				p0Matrix[k][l] = std::sqrt(p0Matrix[k][k] * p0Matrix[l][l]);
				gammaMatrix[k][l] = std::sqrt(
						gammaMatrix[k][k] * gammaMatrix[l][l]);
			}
}

schemi::scalar schemi::mixtureStiffened::Rv() const noexcept
{
	return R;
}

const std::valarray<schemi::scalar>& schemi::mixtureStiffened::Mv() const noexcept
{
	return M;
}

const std::valarray<schemi::scalar>& schemi::mixtureStiffened::Cvv() const noexcept
{
	return CvArr;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::Cv(
		const std::vector<const std::valarray<scalar>*> & concentrations) const noexcept
{
	std::valarray<scalar> CvMixtureOutput(0., concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixtureOutput += X[k] * CvArr[k];

	return CvMixtureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::pFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const noexcept
{
	std::valarray<scalar> pressureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> p0Mixture(0., pressureOutput.size());
	std::valarray<scalar> gammaMixture(0., pressureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];
			gammaMixture += X[k] * X[l] * gammaMatrix[k][l];
		}

	for (std::size_t i = 0; i < pressureOutput.size(); ++i)
		pressureOutput[i] = stiffenedFluid::pFromUv(gammaMixture[i] - 1.0,
				Uv[i], p0Mixture[i], gammaMixture[i]);

	return pressureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::UvFromp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & p) const noexcept
{
	std::valarray<scalar> UvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> p0Mixture(0., UvOutput.size());
	std::valarray<scalar> gammaMixture(0., UvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];
			gammaMixture += X[k] * X[l] * gammaMatrix[k][l];
		}

	for (std::size_t i = 0; i < UvOutput.size(); ++i)
		UvOutput[i] = stiffenedFluid::UvFromp(gammaMixture[i] - 1.0, p[i],
				p0Mixture[i], gammaMixture[i]);

	return UvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::pcFromT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> pressureConcentrationRatioOutput(
			concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> p0Mixture(0.,
			pressureConcentrationRatioOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	for (std::size_t i = 0; i < pressureConcentrationRatioOutput.size(); ++i)
		pressureConcentrationRatioOutput[i] = stiffenedFluid::pcFromT(R, T[i],
				p0Mixture[i], (*concentrations[0])[i]);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::pcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> pressureConcentrationRatioOutput(
			concentration.size());

	const auto & p0k = p0Matrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < pressureConcentrationRatioOutput.size(); ++i)
		pressureConcentrationRatioOutput[i] = stiffenedFluid::pcFromT(R, T[i],
				p0k, concentration[i]);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::UvcFromT(
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

	std::valarray<scalar> p0Mixture(0.,
			InternalEnergyConcentrationRatioOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] = stiffenedFluid::UvcFromT(
				CvMixture[i], T[i], p0Mixture[i], (*concentrations[0])[i]);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::UvcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> InternalEnergyConcentrationRatioOutput(
			concentration.size());

	const auto & Cvk = CvArr[componentIndex];

	const auto & p0k = p0Matrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] = stiffenedFluid::UvcFromT(
				Cvk, T[i], p0k, concentration[i]);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::TFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const noexcept
{
	std::valarray<scalar> TemperatureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., TemperatureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> p0Mixture(0., TemperatureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	for (std::size_t i = 0; i < TemperatureOutput.size(); ++i)
		TemperatureOutput[i] = stiffenedFluid::TFromUv(CvMixture[i], Uv[i],
				(*concentrations[0])[i], p0Mixture[i]);

	return TemperatureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::dpdrho(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> dpdrhoOutput(concentrations[0]->size());

	for (std::size_t i = 0; i < dpdrhoOutput.size(); ++i)
		dpdrhoOutput[i] = stiffenedFluid::dpdrho();

	return dpdrhoOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::dpdUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> dpdUvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> gammaMixture(0., dpdUvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			gammaMixture += X[k] * X[l] * gammaMatrix[k][l];

	for (std::size_t i = 0; i < dpdUvOutput.size(); ++i)
		dpdUvOutput[i] = stiffenedFluid::dpdUv(gammaMixture[i] - 1.0);

	return dpdUvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::nonIdeality(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> nonIdealityOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> p0Mixture(0., nonIdealityOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	for (std::size_t i = 0; i < nonIdealityOutput.size(); ++i)
		nonIdealityOutput[i] = stiffenedFluid::nonIdeality(p0Mixture[i]);

	return nonIdealityOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::sqSonicSpeed(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & density, const std::valarray<scalar> & Uv,
		const std::valarray<scalar> & pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::Cp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> CpMixtureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., concentrations[0]->size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> gammaMixture(0., CpMixtureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			gammaMixture += X[k] * X[l] * gammaMatrix[k][l];

	for (std::size_t i = 0; i < CpMixtureOutput.size(); ++i)
		CpMixtureOutput[i] = stiffenedFluid::Cp(CvMixture[i], gammaMixture[i]);

	return CpMixtureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::Cpk(
		const std::valarray<scalar>&, const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> CpkOutput(T.size());

	for (std::size_t i = 0; i < CpkOutput.size(); ++i)
		CpkOutput[i] = stiffenedFluid::Cp(CvArr[componentIndex],
				gammaMatrix[componentIndex][componentIndex]);

	return CpkOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::hT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::hkT(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::Fv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> FvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	std::valarray<scalar> p0Mixture(0., FvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = stiffenedFluid::Fmx(hPlanck, molecMass, kB, T[i],
				rearX[i], molecVolumeMixture[i], p0Mixture[i]);

	return FvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::Fvk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> FvOutput(concentration.size());

	const auto & Mk = M[componentIndex];

	const auto & p0k = p0Matrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = stiffenedFluid::Fv(concentration[i], T[i], R, Mk, p0k,
				hPlanck);

	return FvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::Sv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> SvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = stiffenedFluid::Smx(hPlanck, molecMass, kB, T[i],
				rearX[i], molecVolumeMixture[i]);

	return SvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureStiffened::Svk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> SvOutput(concentration.size());

	const auto & Mk = M[componentIndex];

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = stiffenedFluid::Sv(concentration[i], T[i], R, Mk,
				hPlanck);

	return SvOutput;
}

schemi::scalar schemi::mixtureStiffened::Cv(
		const std::valarray<scalar> & concentrations) const noexcept
{
	scalar CvMixtureOutput { 0 };

	const auto X = calcMolarFrac(concentrations);

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixtureOutput += X[k] * CvArr[k];

	return CvMixtureOutput;
}

schemi::scalar schemi::mixtureStiffened::pFromUv(
		const std::valarray<scalar> & concentrations,
		const scalar Uv) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar p0Mixture { 0 };
	scalar gammaMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];
			gammaMixture += X[k] * X[l] * gammaMatrix[k][l];
		}

	return stiffenedFluid::pFromUv(gammaMixture - 1.0, Uv, p0Mixture,
			gammaMixture);
}

schemi::scalar schemi::mixtureStiffened::UvFromp(
		const std::valarray<scalar> & concentrations,
		const scalar p) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar p0Mixture { 0 };
	scalar gammaMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];
			gammaMixture += X[k] * X[l] * gammaMatrix[k][l];
		}

	return stiffenedFluid::UvFromp(gammaMixture - 1.0, p, p0Mixture,
			gammaMixture);
}

schemi::scalar schemi::mixtureStiffened::pcFromT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar p0Mixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	return stiffenedFluid::pcFromT(R, T, p0Mixture, concentrations[0]);
}

schemi::scalar schemi::mixtureStiffened::pcFromTk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & p0k = p0Matrix[componentIndex][componentIndex];

	return stiffenedFluid::pcFromT(R, T, p0k, concentration);
}

schemi::scalar schemi::mixtureStiffened::UvcFromT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar p0Mixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	return stiffenedFluid::UvcFromT(CvMixture, T, p0Mixture, concentrations[0]);
}

schemi::scalar schemi::mixtureStiffened::UvcFromTk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & Cvk = CvArr[componentIndex];

	const auto & p0k = p0Matrix[componentIndex][componentIndex];

	return stiffenedFluid::UvcFromT(Cvk, T, p0k, concentration);
}

schemi::scalar schemi::mixtureStiffened::TFromUv(
		const std::valarray<scalar> & concentrations,
		const scalar Uv) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar p0Mixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	return stiffenedFluid::TFromUv(CvMixture, Uv, concentrations[0], p0Mixture);
}

schemi::scalar schemi::mixtureStiffened::dpdrho(const std::valarray<scalar>&,
		const scalar) const noexcept
{
	return stiffenedFluid::dpdrho();
}

schemi::scalar schemi::mixtureStiffened::dpdUv(
		const std::valarray<scalar> & concentrations,
		const scalar) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar gammaMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			gammaMixture += X[k] * X[l] * gammaMatrix[k][l];

	return stiffenedFluid::dpdUv(gammaMixture - 1.0);
}

schemi::scalar schemi::mixtureStiffened::nonIdeality(
		const std::valarray<scalar> & concentrations,
		const scalar) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar p0Mixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	return stiffenedFluid::nonIdeality(p0Mixture);
}

schemi::scalar schemi::mixtureStiffened::sqSonicSpeed(
		const std::valarray<scalar> & concentrations, const scalar density,
		const scalar Uv, const scalar pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

schemi::scalar schemi::mixtureStiffened::Cp(
		const std::valarray<scalar> & concentrations,
		const scalar) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar gammaMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			gammaMixture += X[k] * X[l] * gammaMatrix[k][l];

	return stiffenedFluid::Cp(CvMixture, gammaMixture);
}

schemi::scalar schemi::mixtureStiffened::Cpk(const scalar, const scalar,
		const std::size_t componentIndex) const noexcept
{
	return stiffenedFluid::Cp(CvArr[componentIndex],
			gammaMatrix[componentIndex][componentIndex]);
}

schemi::scalar schemi::mixtureStiffened::hT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

schemi::scalar schemi::mixtureStiffened::hkT(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

schemi::scalar schemi::mixtureStiffened::Fv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	scalar p0Mixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			p0Mixture += X[k] * X[l] * p0Matrix[k][l];

	return stiffenedFluid::Fmx(hPlanck, molecMass, kB, T, X, molecVolumeMixture,
			p0Mixture) * molecConc;
}

schemi::scalar schemi::mixtureStiffened::Fvk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & p0k { p0Matrix[componentIndex][componentIndex] };

	return stiffenedFluid::Fv(concentration, T, R, M[componentIndex], p0k,
			hPlanck);
}

schemi::scalar schemi::mixtureStiffened::Sv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	return stiffenedFluid::Smx(hPlanck, molecMass, kB, T, X, molecVolumeMixture)
			* molecConc;
}

schemi::scalar schemi::mixtureStiffened::Svk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	return stiffenedFluid::Sv(concentration, T, R, M[componentIndex], hPlanck);
}
