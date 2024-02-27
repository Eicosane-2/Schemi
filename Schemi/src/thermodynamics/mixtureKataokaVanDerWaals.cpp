/*
 * mixtureKataokaVanDerWaals.cpp
 *
 *  Created on: 2024/02/23
 *      Author: Maxim Boldyrev
 */

#include "mixtureKataokaVanDerWaals.hpp"

#include "intExpPow.hpp"
#include "globalConstants.hpp"

schemi::mixtureKataokaVanDerWaals::mixtureKataokaVanDerWaals() noexcept :
		abstractMixtureThermodynamics(0, 0), M(0), CvArr(0), molecMass(0), Tcrit(
				0), Pcrit(0), Vcrit(0), V0Matrix(0), epsMatrix(0), bMatrix(0), V0MatrixMolec(
				0), epsMatrixMolec(0), bMatrixMolec(0)
{
}

schemi::mixtureKataokaVanDerWaals::mixtureKataokaVanDerWaals(const scalar Rin,
		const scalar hPin, const std::valarray<scalar> & Min,
		const std::valarray<scalar> & Cvin,
		const std::valarray<scalar> & Tcritin,
		const std::valarray<scalar> & Pcritin,
		const std::valarray<scalar> & epsilonLJ, /*J/mole*/
		const std::valarray<scalar> & sigmaLJ, /*m*/
		const std::pair<bool, scalar> & bCalcType) noexcept :
		abstractMixtureThermodynamics(Rin, hPin), M(Min), CvArr(Cvin), molecMass(
				M / NAvogardro), Tcrit(Tcritin), Pcrit(Pcritin), Vcrit(
				Min.size()), V0Matrix(std::valarray<scalar>(Min.size()),
				Min.size()), epsMatrix(std::valarray<scalar>(Min.size()),
				Min.size()), bMatrix(std::valarray<scalar>(Min.size()),
				Min.size()), V0MatrixMolec(std::valarray<scalar>(Min.size()),
				Min.size()), epsMatrixMolec(std::valarray<scalar>(Min.size()),
				Min.size()), bMatrixMolec(std::valarray<scalar>(Min.size()),
				Min.size())
{
	const std::size_t numberOfComponents { Min.size() };

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		bMatrix[k][k] = [this, &bCalcType, &sigmaLJ](const std::size_t comp)
		{
			if (bCalcType.first)
				return bCalcType.second * R * Tcrit[comp] / (8 * Pcrit[comp]);
			else
				return bCalcType.second * pow<scalar, 3>(sigmaLJ[comp])
						* NAvogardro;
		}(k);
		V0Matrix[k][k] = pow<scalar, 3>(sigmaLJ[k]) / std::sqrt(2) * NAvogardro;
		epsMatrix[k][k] = epsilonLJ[k];

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
				V0Matrix[k][l] = std::pow(
						0.5
								* (std::cbrt(V0Matrix[k][k])
										+ std::cbrt(V0Matrix[l][l])), 3);
				epsMatrix[k][l] = std::sqrt(epsMatrix[k][k] * epsMatrix[l][l]);
			}

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		bMatrixMolec[k] = bMatrix[k] / NAvogardro;
		V0MatrixMolec[k] = V0Matrix[k] / NAvogardro;
		epsMatrixMolec[k] = epsMatrix[k] / NAvogardro;
	}
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::Rv() const noexcept
{
	return R;
}

const std::valarray<schemi::scalar>& schemi::mixtureKataokaVanDerWaals::Mv() const noexcept
{
	return M;
}

const std::valarray<schemi::scalar>& schemi::mixtureKataokaVanDerWaals::Cvv() const noexcept
{
	return CvArr;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::Cv(
		const std::vector<const std::valarray<scalar>*> & concentrations) const noexcept
{
	std::valarray<scalar> CvMixOutput(0., concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixOutput += X[k] * CvArr[k];

	return CvMixOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::pFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const
{
	std::valarray<scalar> pressureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., pressureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> V0Mixture(0., pressureOutput.size()), epsMixture(0.,
			pressureOutput.size()), bMixture(0., pressureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < pressureOutput.size(); ++i)
		pressureOutput[i] = KataokaVanDerWaalsFluid::pFromUv(
				(*concentrations[0])[i], CvMixture[i], Uv[i], V0Mixture[i],
				epsMixture[i], R, bMixture[i]);

	return pressureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::UvFromp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & p) const
{
	std::valarray<scalar> UvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., UvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> V0Mixture(0., UvOutput.size()), epsMixture(0.,
			UvOutput.size()), bMixture(0., UvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < UvOutput.size(); ++i)
		UvOutput[i] = KataokaVanDerWaalsFluid::UvFromp((*concentrations[0])[i],
				CvMixture[i], p[i], V0Mixture[i], epsMixture[i], R,
				bMixture[i]);

	return UvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::pcFromT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> pressureConcentrationRatioOutput(
			concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> V0Mixture(0.,
			pressureConcentrationRatioOutput.size()), epsMixture(0.,
			pressureConcentrationRatioOutput.size()), bMixture(0.,
			pressureConcentrationRatioOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < pressureConcentrationRatioOutput.size(); ++i)
		pressureConcentrationRatioOutput[i] = KataokaVanDerWaalsFluid::pcFromT(
				(*concentrations[0])[i], T[i], V0Mixture[i], epsMixture[i], R,
				bMixture[i]);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::pcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> pressureConcentrationRatioOutput(
			concentration.size());

	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < pressureConcentrationRatioOutput.size(); ++i)
		pressureConcentrationRatioOutput[i] = KataokaVanDerWaalsFluid::pcFromT(
				concentration[i], T[i], V0k, epsk, R, bk);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::UvcFromT(
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

	std::valarray<scalar> V0Mixture(0.,
			InternalEnergyConcentrationRatioOutput.size()), epsMixture(0.,
			InternalEnergyConcentrationRatioOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
		}

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] =
				KataokaVanDerWaalsFluid::UvcFromT((*concentrations[0])[i],
						CvMixture[i], T[i], V0Mixture[i], epsMixture[i], R);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::UvcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> InternalEnergyConcentrationRatioOutput(
			concentration.size());

	const auto & Cvk = CvArr[componentIndex];

	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] =
				KataokaVanDerWaalsFluid::UvcFromT(concentration[i], Cvk, T[i],
						V0k, epsk, R);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::TFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const
{
	std::valarray<scalar> TemperatureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., TemperatureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> V0Mixture(0., TemperatureOutput.size()), epsMixture(
			0., TemperatureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
		}

	for (std::size_t i = 0; i < TemperatureOutput.size(); ++i)
		TemperatureOutput[i] = KataokaVanDerWaalsFluid::TFromUv(
				(*concentrations[0])[i], CvMixture[i], Uv[i], V0Mixture[i],
				epsMixture[i], R);

	return TemperatureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::dpdrho(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const
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

	std::valarray<scalar> V0Mixture(0., dpdrhoOutput.size()), epsMixture(0.,
			dpdrhoOutput.size()), bMixture(0., dpdrhoOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < dpdrhoOutput.size(); ++i)
		dpdrhoOutput[i] = KataokaVanDerWaalsFluid::dpdrho(MMixture[i],
				(*concentrations[0])[i], CvMixture[i], Uv[i], V0Mixture[i],
				epsMixture[i], R, bMixture[i]);

	return dpdrhoOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::dpdUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const
{
	std::valarray<scalar> dpdUvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., dpdUvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> V0Mixture(0., dpdUvOutput.size()), epsMixture(0.,
			dpdUvOutput.size()), bMixture(0., dpdUvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < dpdUvOutput.size(); ++i)
		dpdUvOutput[i] = KataokaVanDerWaalsFluid::dpdUv((*concentrations[0])[i],
				CvMixture[i], Uv[i], V0Mixture[i], epsMixture[i], R,
				bMixture[i]);

	return dpdUvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::nonIdeality(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> nonIdealityOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> V0Mixture(0., nonIdealityOutput.size()), epsMixture(
			0., nonIdealityOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
		}

	for (std::size_t i = 0; i < nonIdealityOutput.size(); ++i)
		nonIdealityOutput[i] = KataokaVanDerWaalsFluid::nonIdeality(
				(*concentrations[0])[i], T[i], V0Mixture[i], epsMixture[i], R);

	return nonIdealityOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::sqSonicSpeed(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & density, const std::valarray<scalar> & Uv,
		const std::valarray<scalar> & pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::Cp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> CpMixOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., concentrations[0]->size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> V0Mixture(0., CpMixOutput.size()), epsMixture(0.,
			CpMixOutput.size()), bMixture(0., CpMixOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < CpMixOutput.size(); ++i)
		CpMixOutput[i] = KataokaVanDerWaalsFluid::Cp((*concentrations[0])[i],
				T[i], V0Mixture[i], epsMixture[i], bMixture[i], R,
				CvMixture[i]);

	return CpMixOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::Cpk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> CpkOutput(T.size());

	const auto & Cvk = CvArr[componentIndex];

	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < CpkOutput.size(); ++i)
		CpkOutput[i] = KataokaVanDerWaalsFluid::Cp(concentration[i], T[i], V0k,
				epsk, bk, R, Cvk);

	return CpkOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::hT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::hkT(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::Fv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> FvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	std::valarray<scalar> V0MixtureMolec(0., FvOutput.size()), epsMixtureMolec(
			0., FvOutput.size()), bMixtureMolec(0., FvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0MixtureMolec += X[k] * X[l] * V0MatrixMolec[k][l];
			epsMixtureMolec += X[k] * X[l] * epsMatrixMolec[k][l];
			bMixtureMolec += X[k] * X[l] * bMatrixMolec[k][l];
		}

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = KataokaVanDerWaalsFluid::Fmx(hPlanck, molecMass, kB, T[i],
				rearX[i], molecVolumeMixture[i], V0MixtureMolec[i],
				epsMixtureMolec[i], bMixtureMolec[i]);

	return FvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::Fvk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> FvOutput(concentration.size());

	const auto & Mk = M[componentIndex];

	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = KataokaVanDerWaalsFluid::Fv(concentration[i], T[i], R, Mk,
				V0k, epsk, bk, hPlanck);

	return FvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::Sv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> SvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	std::valarray<scalar> V0MixtureMolec(0., SvOutput.size()), epsMixtureMolec(
			0., SvOutput.size()), bMixtureMolec(0., SvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0MixtureMolec += X[k] * X[l] * V0MatrixMolec[k][l];
			epsMixtureMolec += X[k] * X[l] * epsMatrixMolec[k][l];
			bMixtureMolec += X[k] * X[l] * bMatrixMolec[k][l];
		}

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = KataokaVanDerWaalsFluid::Smx(hPlanck, molecMass, kB, T[i],
				rearX[i], molecVolumeMixture[i], V0MixtureMolec[i],
				epsMixtureMolec[i], bMixtureMolec[i]);

	return SvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureKataokaVanDerWaals::Svk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> SvOutput(concentration.size());

	const auto & Mk = M[componentIndex];

	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = KataokaVanDerWaalsFluid::Sv(concentration[i], T[i], R, Mk,
				V0k, epsk, bk, hPlanck);

	return SvOutput;
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::Cv(
		const std::valarray<scalar> & concentrations) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixOutput { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixOutput += X[k] * CvArr[k];

	return CvMixOutput;
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::pFromUv(
		const std::valarray<scalar> & concentrations, const scalar Uv) const
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar V0Mixture { 0. }, epsMixture { 0. }, bMixture { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return KataokaVanDerWaalsFluid::pFromUv(concentrations[0], CvMixture, Uv,
			V0Mixture, epsMixture, R, bMixture);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::UvFromp(
		const std::valarray<scalar> & concentrations, const scalar p) const
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar V0Mixture { 0. }, epsMixture { 0. }, bMixture { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return KataokaVanDerWaalsFluid::UvFromp(concentrations[0], CvMixture, p,
			V0Mixture, epsMixture, R, bMixture);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::pcFromT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar V0Mixture { 0. }, epsMixture { 0. }, bMixture { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return KataokaVanDerWaalsFluid::pcFromT(concentrations[0], T, V0Mixture,
			epsMixture, R, bMixture);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::pcFromTk(
		const scalar concentration, const scalar T,
		const std::size_t componentIndex) const noexcept
{
	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	return KataokaVanDerWaalsFluid::pcFromT(concentration, T, V0k, epsk, R, bk);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::UvcFromT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar V0Mixture { 0. }, epsMixture { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
		}

	return KataokaVanDerWaalsFluid::UvcFromT(concentrations[0], CvMixture, T,
			V0Mixture, epsMixture, R);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::UvcFromTk(
		const scalar concentration, const scalar T,
		const std::size_t componentIndex) const noexcept
{
	const auto & Cvk = CvArr[componentIndex];

	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];

	return KataokaVanDerWaalsFluid::UvcFromT(concentration, Cvk, T, V0k, epsk,
			R);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::TFromUv(
		const std::valarray<scalar> & concentrations, const scalar Uv) const
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar V0Mixture { 0. }, epsMixture { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
		}

	return KataokaVanDerWaalsFluid::TFromUv(concentrations[0], CvMixture, Uv,
			V0Mixture, epsMixture, R);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::dpdrho(
		const std::valarray<scalar> & concentrations, const scalar Uv) const
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 }, MMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
	{
		CvMixture += X[k] * CvArr[k];
		MMixture += X[k] * M[k];
	}

	scalar V0Mixture { 0. }, epsMixture { 0. }, bMixture { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return KataokaVanDerWaalsFluid::dpdrho(MMixture, concentrations[0],
			CvMixture, Uv, V0Mixture, epsMixture, R, bMixture);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::dpdUv(
		const std::valarray<scalar> & concentrations, const scalar Uv) const
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar V0Mixture { 0. }, epsMixture { 0. }, bMixture { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return KataokaVanDerWaalsFluid::dpdUv(concentrations[0], CvMixture, Uv,
			V0Mixture, epsMixture, R, bMixture);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::nonIdeality(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar V0Mixture { 0. }, epsMixture { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
		}

	return KataokaVanDerWaalsFluid::nonIdeality(concentrations[0], T, V0Mixture,
			epsMixture, R);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::sqSonicSpeed(
		const std::valarray<scalar> & concentrations, const scalar density,
		const scalar Uv, const scalar pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::Cp(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixture { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	scalar V0Mixture { 0. }, epsMixture { 0. }, bMixture { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0Mixture += X[k] * X[l] * V0Matrix[k][l];
			epsMixture += X[k] * X[l] * epsMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	return KataokaVanDerWaalsFluid::Cp(concentrations[0], T, V0Mixture,
			epsMixture, bMixture, R, CvMixture);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::Cpk(
		const scalar concentration, const scalar T,
		const std::size_t componentIndex) const noexcept
{
	const auto & Cvk = CvArr[componentIndex];

	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	return KataokaVanDerWaalsFluid::Cp(concentration, T, V0k, epsk, bk, R, Cvk);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::hT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::hkT(
		const scalar concentration, const scalar T,
		const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::Fv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	scalar V0MixtureMolec { 0. }, epsMixtureMolec { 0. }, bMixtureMolec { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0MixtureMolec += X[k] * X[l] * V0MatrixMolec[k][l];
			epsMixtureMolec += X[k] * X[l] * epsMatrixMolec[k][l];
			bMixtureMolec += X[k] * X[l] * bMatrixMolec[k][l];
		}

	return KataokaVanDerWaalsFluid::Fmx(hPlanck, molecMass, kB, T, X,
			molecVolumeMixture, V0MixtureMolec, epsMixtureMolec, bMixtureMolec)
			* molecConc;
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::Fvk(
		const scalar concentration, const scalar T,
		const std::size_t componentIndex) const noexcept
{
	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	return KataokaVanDerWaalsFluid::Fv(concentration, T, R, M[componentIndex],
			V0k, epsk, bk, hPlanck);
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::Sv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	scalar V0MixtureMolec { 0. }, epsMixtureMolec { 0. }, bMixtureMolec { 0. };

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			V0MixtureMolec += X[k] * X[l] * V0MatrixMolec[k][l];
			epsMixtureMolec += X[k] * X[l] * epsMatrixMolec[k][l];
			bMixtureMolec += X[k] * X[l] * bMatrixMolec[k][l];
		}

	return KataokaVanDerWaalsFluid::Smx(hPlanck, molecMass, kB, T, X,
			molecVolumeMixture, V0MixtureMolec, epsMixtureMolec, bMixtureMolec)
			* molecConc;
}

schemi::scalar schemi::mixtureKataokaVanDerWaals::Svk(
		const scalar concentration, const scalar T,
		const std::size_t componentIndex) const noexcept
{
	const auto & V0k = V0Matrix[componentIndex][componentIndex];
	const auto & epsk = epsMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	return KataokaVanDerWaalsFluid::Sv(concentration, T, R, M[componentIndex],
			V0k, epsk, bk, hPlanck);
}
