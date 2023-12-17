/*
 * mixtureVanDerWaals.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "mixtureVanDerWaals.hpp"

#include "globalConstants.hpp"

schemi::mixtureVanDerWaals::mixtureVanDerWaals() noexcept :
		abstractMixtureThermodynamics(0, 0), M(0), CvArr(0), molecMass(0), Tcrit(
				0), Pcrit(0), Vcrit(0), aMatrix(0), bMatrix(0)
{
}

schemi::mixtureVanDerWaals::mixtureVanDerWaals(const scalar Rin,
		const scalar hPin, const std::valarray<scalar> & Min,
		const std::valarray<scalar> & Cvin,
		const std::valarray<scalar> & Tcritin,
		const std::valarray<scalar> & Pcritin) noexcept :
		abstractMixtureThermodynamics(Rin, hPin), M(Min), CvArr(Cvin), molecMass(
				M / NAvogardro), Tcrit(Tcritin), Pcrit(Pcritin), Vcrit(
				Min.size()), aMatrix(std::valarray<scalar>(Min.size()),
				Min.size()), bMatrix(std::valarray<scalar>(Min.size()),
				Min.size())
{
	const std::size_t numberOfComponents { Min.size() };

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		bMatrix[k][k] = Tcrit[k] * R / (8 * Pcrit[k]);
		aMatrix[k][k] = Pcrit[k] * bMatrix[k][k] * bMatrix[k][k] * 27;

		Vcrit[k] = 3 * bMatrix[k][k];
	}

	for (std::size_t k = 0; k < numberOfComponents; ++k)
		for (std::size_t l = 0; l < numberOfComponents; ++l)
			if (k != l)
			{
				bMatrix[k][l] = std::pow(
						0.5
								* (std::cbrt(bMatrix[k][k])
										+ std::cbrt(bMatrix[l][l])), 3);
				aMatrix[k][l] = std::sqrt(aMatrix[k][k] * aMatrix[l][l]);
			}
}

schemi::scalar schemi::mixtureVanDerWaals::Rv() const noexcept
{
	return R;
}

const std::valarray<schemi::scalar>& schemi::mixtureVanDerWaals::Mv() const noexcept
{
	return M;
}

const std::valarray<schemi::scalar>& schemi::mixtureVanDerWaals::Cvv() const noexcept
{
	return CvArr;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::Cv(
		const std::vector<const std::valarray<scalar>*> & concentrations) const noexcept
{
	std::valarray<scalar> CvMixtureOutput(0., concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixtureOutput += X[k] * CvArr[k];

	return CvMixtureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::pFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const noexcept
{
	std::valarray<scalar> pressureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., pressureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> aMixture(0., pressureOutput.size());
	std::valarray<scalar> bMixture(0., pressureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < pressureOutput.size(); ++i)
		pressureOutput[i] = vanDerWaalsFluid::pFromUv(R / CvMixture[i], Uv[i],
				(*concentrations[0])[i], aMixture[i], bMixture[i]);

	return pressureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::UvFromp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & p) const noexcept
{
	std::valarray<scalar> UvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., UvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> aMixture(0., UvOutput.size());
	std::valarray<scalar> bMixture(0., UvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < UvOutput.size(); ++i)
		UvOutput[i] = vanDerWaalsFluid::UvFromp(R / CvMixture[i], p[i],
				(*concentrations[0])[i], aMixture[i], bMixture[i]);

	return UvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::pcFromT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> pressureConcentrationRatioOutput(
			concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> aMixture(0., pressureConcentrationRatioOutput.size());
	std::valarray<scalar> bMixture(0., pressureConcentrationRatioOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < pressureConcentrationRatioOutput.size(); ++i)
		pressureConcentrationRatioOutput[i] = vanDerWaalsFluid::pcFromT(R, T[i],
				(*concentrations[0])[i], aMixture[i], bMixture[i]);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::pcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> pressureConcentrationRatioOutput(
			concentration.size());

	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < pressureConcentrationRatioOutput.size(); ++i)
		pressureConcentrationRatioOutput[i] = vanDerWaalsFluid::pcFromT(R, T[i],
				concentration[i], ak, bk);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::UvcFromT(
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

	std::valarray<scalar> aMixture(0.,
			InternalEnergyConcentrationRatioOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			aMixture += X[k] * X[l] * aMatrix[k][l];

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] = vanDerWaalsFluid::UvcFromT(
				CvMixture[i], T[i], (*concentrations[0])[i], aMixture[i]);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::UvcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> InternalEnergyConcentrationRatioOutput(
			concentration.size());

	const auto & Cvk = CvArr[componentIndex];

	const auto & ak = aMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] = vanDerWaalsFluid::UvcFromT(
				Cvk, T[i], concentration[i], ak);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::TFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const noexcept
{
	std::valarray<scalar> TemperatureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., TemperatureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> aMixture(0., TemperatureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			aMixture += X[k] * X[l] * aMatrix[k][l];

	for (std::size_t i = 0; i < TemperatureOutput.size(); ++i)
		TemperatureOutput[i] = vanDerWaalsFluid::TFromUv(CvMixture[i], Uv[i],
				(*concentrations[0])[i], aMixture[i]);

	return TemperatureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::dpdrho(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const noexcept
{
	std::valarray<scalar> dpdrhoOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., dpdrhoOutput.size());
	std::valarray<scalar> MMixture(0., dpdrhoOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
	{
		CvMixture += X[k] * CvArr[k];
		MMixture += X[k] * M[k];
	}

	std::valarray<scalar> aMixture(0., dpdrhoOutput.size());
	std::valarray<scalar> bMixture(0., dpdrhoOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < dpdrhoOutput.size(); ++i)
		dpdrhoOutput[i] = vanDerWaalsFluid::dpdrho(R / CvMixture[i], Uv[i],
				(*concentrations[0])[i], MMixture[i], aMixture[i], bMixture[i]);

	return dpdrhoOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::dpdUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> dpdUvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., dpdUvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> bMixture(0., dpdUvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			bMixture += X[k] * X[l] * bMatrix[k][l];

	for (std::size_t i = 0; i < dpdUvOutput.size(); ++i)
		dpdUvOutput[i] = vanDerWaalsFluid::dpdUv(R / CvMixture[i],
				(*concentrations[0])[i], bMixture[i]);

	return dpdUvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::nonIdeality(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar>&) const noexcept
{
	std::valarray<scalar> nonIdealityOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> aMixture(0., nonIdealityOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			aMixture += X[k] * X[l] * aMatrix[k][l];

	for (std::size_t i = 0; i < nonIdealityOutput.size(); ++i)
		nonIdealityOutput[i] = vanDerWaalsFluid::nonIdeality(
				(*concentrations[0])[i], aMixture[i]);

	return nonIdealityOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::sqSonicSpeed(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & density, const std::valarray<scalar> & Uv,
		const std::valarray<scalar> & pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::Cp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> CpMixtureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., concentrations[0]->size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> aMixture(0., CpMixtureOutput.size());
	std::valarray<scalar> bMixture(0., CpMixtureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < CpMixtureOutput.size(); ++i)
		CpMixtureOutput[i] = vanDerWaalsFluid::Cp(CvMixture[i], R, aMixture[i],
				bMixture[i], (*concentrations[0])[i], T[i]);

	return CpMixtureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::Cpk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> CpkOutput(T.size());

	for (std::size_t i = 0; i < CpkOutput.size(); ++i)
		CpkOutput[i] = vanDerWaalsFluid::Cp(CvArr[componentIndex], R,
				aMatrix[componentIndex][componentIndex],
				bMatrix[componentIndex][componentIndex], concentration[i],
				T[i]);

	return CpkOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::hT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::hkT(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::Fv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> FvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	std::valarray<scalar> aMixture(0., FvOutput.size());
	std::valarray<scalar> bMixture(0., FvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = vanDerWaalsFluid::Fmx(hPlanck, molecMass, kB, T[i],
				rearX[i], molecVolumeMixture[i],
				aMixture[i] / (NAvogardro * NAvogardro),
				bMixture[i] / NAvogardro);

	return FvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::Fvk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> FvOutput(concentration.size());

	const auto & Mk = M[componentIndex];

	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = vanDerWaalsFluid::Fv(concentration[i], T[i], R, Mk, ak,
				bk, hPlanck);

	return FvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::Sv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> SvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	std::valarray<scalar> bMixture(0., SvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			bMixture += X[k] * X[l] * bMatrix[k][l];

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = vanDerWaalsFluid::Smx(hPlanck, molecMass, kB, T[i],
				rearX[i], molecVolumeMixture[i], bMixture[i] / NAvogardro);

	return SvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureVanDerWaals::Svk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> SvOutput(concentration.size());

	const auto & Mk = M[componentIndex];

	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = vanDerWaalsFluid::Sv(concentration[i], T[i], R, Mk, bk,
				hPlanck);

	return SvOutput;
}

schemi::scalar schemi::mixtureVanDerWaals::Cv(
		const std::valarray<scalar> & concentrations) const noexcept
{
	scalar CvMixtureOutput { 0 };

	const auto X = calcMolarFrac(concentrations);

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixtureOutput += X[k] * CvArr[k];

	return CvMixtureOutput;
}

schemi::scalar schemi::mixtureVanDerWaals::pFromUv(
		const std::valarray<scalar> & concentrations,
		const scalar Uv) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar aMixture { 0 };
	scalar bMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return vanDerWaalsFluid::pFromUv(R / CvMixture, Uv, concentrations[0],
			aMixture, bMixture);
}

schemi::scalar schemi::mixtureVanDerWaals::UvFromp(
		const std::valarray<scalar> & concentrations,
		const scalar p) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar aMixture { 0 };
	scalar bMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return vanDerWaalsFluid::UvFromp(R / CvMixture, p, concentrations[0],
			aMixture, bMixture);
}

schemi::scalar schemi::mixtureVanDerWaals::pcFromT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar aMixture { 0 };
	scalar bMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return vanDerWaalsFluid::pcFromT(R, T, concentrations[0], aMixture,
			bMixture);
}

schemi::scalar schemi::mixtureVanDerWaals::pcFromTk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	return vanDerWaalsFluid::pcFromT(R, T, concentration, ak, bk);
}

schemi::scalar schemi::mixtureVanDerWaals::UvcFromT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar aMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			aMixture += X[k] * X[l] * aMatrix[k][l];

	return vanDerWaalsFluid::UvcFromT(CvMixture, T, concentrations[0], aMixture);
}

schemi::scalar schemi::mixtureVanDerWaals::UvcFromTk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & Cvk = CvArr[componentIndex];

	const auto & ak = aMatrix[componentIndex][componentIndex];

	return vanDerWaalsFluid::UvcFromT(Cvk, T, concentration, ak);
}

schemi::scalar schemi::mixtureVanDerWaals::TFromUv(
		const std::valarray<scalar> & concentrations,
		const scalar Uv) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar aMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			aMixture += X[k] * X[l] * aMatrix[k][l];

	return vanDerWaalsFluid::TFromUv(CvMixture, Uv, concentrations[0], aMixture);
}

schemi::scalar schemi::mixtureVanDerWaals::dpdrho(
		const std::valarray<scalar> & concentrations,
		const scalar Uv) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };
	scalar MMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
	{
		CvMixture += X[k] * CvArr[k];
		MMixture += X[k] * M[k];
	}

	scalar aMixture { 0 };
	scalar bMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return vanDerWaalsFluid::dpdrho(R / CvMixture, Uv, concentrations[0],
			MMixture, aMixture, bMixture);
}

schemi::scalar schemi::mixtureVanDerWaals::dpdUv(
		const std::valarray<scalar> & concentrations,
		const scalar) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar bMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			bMixture += X[k] * X[l] * bMatrix[k][l];

	return vanDerWaalsFluid::dpdUv(R / CvMixture, concentrations[0], bMixture);
}

schemi::scalar schemi::mixtureVanDerWaals::nonIdeality(
		const std::valarray<scalar> & concentrations,
		const scalar) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar aMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			aMixture += X[k] * X[l] * aMatrix[k][l];

	return vanDerWaalsFluid::nonIdeality(concentrations[0], aMixture);
}

schemi::scalar schemi::mixtureVanDerWaals::sqSonicSpeed(
		const std::valarray<scalar> & concentrations, const scalar density,
		const scalar Uv, const scalar pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

schemi::scalar schemi::mixtureVanDerWaals::Cp(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar aMixture { 0 };
	scalar bMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return vanDerWaalsFluid::Cp(CvMixture, R, aMixture, bMixture,
			concentrations[0], T);
}

schemi::scalar schemi::mixtureVanDerWaals::Cpk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	return vanDerWaalsFluid::Cp(CvArr[componentIndex], R,
			aMatrix[componentIndex][componentIndex],
			bMatrix[componentIndex][componentIndex], concentration, T);
}

schemi::scalar schemi::mixtureVanDerWaals::hT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

schemi::scalar schemi::mixtureVanDerWaals::hkT(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

schemi::scalar schemi::mixtureVanDerWaals::Fv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	scalar aMixture { 0 };
	scalar bMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return vanDerWaalsFluid::Fmx(hPlanck, molecMass, kB, T, X,
			molecVolumeMixture, aMixture / (NAvogardro * NAvogardro),
			bMixture / NAvogardro) * molecConc;
}

schemi::scalar schemi::mixtureVanDerWaals::Fvk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & ak { aMatrix[componentIndex][componentIndex] };
	const auto & bk { bMatrix[componentIndex][componentIndex] };

	return vanDerWaalsFluid::Fv(concentration, T, R, M[componentIndex], ak, bk,
			hPlanck);
}

schemi::scalar schemi::mixtureVanDerWaals::Sv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	scalar bMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
			bMixture += X[k] * X[l] * bMatrix[k][l];

	return vanDerWaalsFluid::Smx(hPlanck, molecMass, kB, T, X,
			molecVolumeMixture, bMixture / NAvogardro) * molecConc;
}

schemi::scalar schemi::mixtureVanDerWaals::Svk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & bk { bMatrix[componentIndex][componentIndex] };

	return vanDerWaalsFluid::Sv(concentration, T, R, M[componentIndex], bk,
			hPlanck);
}
