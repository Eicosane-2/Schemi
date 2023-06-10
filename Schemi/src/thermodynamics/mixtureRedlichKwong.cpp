/*
 * mixtureRedlichKwong.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "mixtureRedlichKwong.hpp"

#include "globalConstants.hpp"

schemi::mixtureRedlichKwong::mixtureRedlichKwong() noexcept :
		abstractMixtureThermodynamics(0, 0), M(0), CvArr(0), molecMass(0), Tcrit(
				0), Pcrit(0), Vcrit(0), aMatrix(0), bMatrix(0), aMatrixMolec(0), bMatrixMolec(
				0)
{
}

schemi::mixtureRedlichKwong::mixtureRedlichKwong(const scalar Rin,
		const scalar hPin, const std::valarray<scalar> & Min,
		const std::valarray<scalar> & Cvin,
		const std::valarray<scalar> & Tcritin,
		const std::valarray<scalar> & Pcritin) noexcept :
		abstractMixtureThermodynamics(Rin, hPin), M(Min), CvArr(Cvin), molecMass(
				M / NAvogardro), Tcrit(Tcritin), Pcrit(Pcritin), Vcrit(
				Min.size()), aMatrix(std::valarray<scalar>(Min.size()),
				Min.size()), bMatrix(std::valarray<scalar>(Min.size()),
				Min.size()), aMatrixMolec(std::valarray<scalar>(Min.size()),
				Min.size()), bMatrixMolec(std::valarray<scalar>(Min.size()),
				Min.size())
{
#if ( defined(__GNUG__) ) && ( !defined(__ICC) )
	constexpr scalar chi { std::cbrt(2) - 1. };
	constexpr scalar Sigma_a { 1. / (9. * chi) }, Sigma_b { chi / 3. };
#else
		const scalar chi { std::cbrt(2) - 1. };
		const scalar Sigma_a { 1. / (9. * chi) }, Sigma_b { chi / 3. };
#endif

	const std::size_t numberOfComponents { Min.size() };

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		bMatrix[k][k] = Sigma_b * R * Tcrit[k] / Pcrit[k];
		aMatrix[k][k] = Sigma_a * R * R * Tcrit[k] * Tcrit[k]
				* std::sqrt(Tcrit[k]) / Pcrit[k];

		Vcrit[k] = 1 / 0.259921 * bMatrix[k][k];
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

	for (std::size_t k = 0; k < numberOfComponents; ++k)
	{
		bMatrixMolec[k] = bMatrix[k] / NAvogardro;
		aMatrixMolec[k] = aMatrix[k] / (NAvogardro * NAvogardro);
	}
}

schemi::scalar schemi::mixtureRedlichKwong::Rv() const noexcept
{
	return R;
}

const std::valarray<schemi::scalar>& schemi::mixtureRedlichKwong::Mv() const noexcept
{
	return M;
}

const std::valarray<schemi::scalar>& schemi::mixtureRedlichKwong::Cvv() const noexcept
{
	return CvArr;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::Cv(
		const std::vector<const std::valarray<scalar>*> & concentrations) const noexcept
{
	std::valarray<scalar> CvMixOutput(0., concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixOutput += X[k] * CvArr[k];

	return CvMixOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::pFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const
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
		pressureOutput[i] = RedlichKwongFluid::pFromUv(R, CvMixture[i], Uv[i],
				(*concentrations[0])[i], aMixture[i], bMixture[i]);

	return pressureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::UvFromp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & p) const
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
		UvOutput[i] = RedlichKwongFluid::UvFromp(R, CvMixture[i], p[i],
				(*concentrations[0])[i], aMixture[i], bMixture[i]);

	return UvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::pcFromT(
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
		pressureConcentrationRatioOutput[i] = RedlichKwongFluid::pcFromT(R,
				T[i], (*concentrations[0])[i], aMixture[i], bMixture[i]);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::pcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> pressureConcentrationRatioOutput(
			concentration.size());

	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < pressureConcentrationRatioOutput.size(); ++i)
		pressureConcentrationRatioOutput[i] = RedlichKwongFluid::pcFromT(R,
				T[i], concentration[i], ak, bk);

	return pressureConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::UvcFromT(
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
	std::valarray<scalar> bMixture(0.,
			InternalEnergyConcentrationRatioOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] = RedlichKwongFluid::UvcFromT(
				CvMixture[i], T[i], (*concentrations[0])[i], aMixture[i],
				bMixture[i]);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::UvcFromTk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> InternalEnergyConcentrationRatioOutput(
			concentration.size());

	const auto & Cvk = CvArr[componentIndex];

	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < InternalEnergyConcentrationRatioOutput.size();
			++i)
		InternalEnergyConcentrationRatioOutput[i] = RedlichKwongFluid::UvcFromT(
				Cvk, T[i], concentration[i], ak, bk);

	return InternalEnergyConcentrationRatioOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::TFromUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const
{
	std::valarray<scalar> TemperatureOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., TemperatureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> aMixture(0., TemperatureOutput.size());
	std::valarray<scalar> bMixture(0., TemperatureOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < TemperatureOutput.size(); ++i)
		TemperatureOutput[i] = RedlichKwongFluid::TFromUv(CvMixture[i], Uv[i],
				(*concentrations[0])[i], aMixture[i], bMixture[i]);

	return TemperatureOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::dpdrho(
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

	std::valarray<scalar> aMixture(0., dpdrhoOutput.size());
	std::valarray<scalar> bMixture(0., dpdrhoOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < dpdrhoOutput.size(); ++i)
		dpdrhoOutput[i] = RedlichKwongFluid::dpdrho(CvMixture[i],
				R / CvMixture[i], Uv[i], (*concentrations[0])[i], MMixture[i],
				aMixture[i], bMixture[i]);

	return dpdrhoOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::dpdUv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & Uv) const
{
	std::valarray<scalar> dpdUvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., dpdUvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> aMixture(0., dpdUvOutput.size());
	std::valarray<scalar> bMixture(0., dpdUvOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < dpdUvOutput.size(); ++i)
		dpdUvOutput[i] = RedlichKwongFluid::dpdUv(CvMixture[i],
				R / CvMixture[i], (*concentrations[0])[i], Uv[i], aMixture[i],
				bMixture[i]);

	return dpdUvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::nonIdeality(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> nonIdealityOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> aMixture(0., nonIdealityOutput.size());
	std::valarray<scalar> bMixture(0., nonIdealityOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < nonIdealityOutput.size(); ++i)
		nonIdealityOutput[i] = RedlichKwongFluid::nonIdeality(
				(*concentrations[0])[i], aMixture[i], bMixture[i], T[i]);

	return nonIdealityOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::sqSonicSpeed(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & density, const std::valarray<scalar> & Uv,
		const std::valarray<scalar> & pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::Cp(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> CpMixOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	std::valarray<scalar> CvMixture(0., concentrations[0]->size());

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixture += X[k] * CvArr[k];

	std::valarray<scalar> aMixture(0., CpMixOutput.size());
	std::valarray<scalar> bMixture(0., CpMixOutput.size());

	for (std::size_t k = 0; k < X.size(); ++k)
		for (std::size_t l = 0; l < X.size(); ++l)
		{
			aMixture += X[k] * X[l] * aMatrix[k][l];
			bMixture += X[k] * X[l] * bMatrix[k][l];
		}

	for (std::size_t i = 0; i < CpMixOutput.size(); ++i)
		CpMixOutput[i] = RedlichKwongFluid::Cp(CvMixture[i], R, aMixture[i],
				bMixture[i], (*concentrations[0])[i], T[i]);

	return CpMixOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::Cpk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> CpkOutput(T.size());

	for (std::size_t i = 0; i < CpkOutput.size(); ++i)
		CpkOutput[i] = RedlichKwongFluid::Cp(CvArr[componentIndex], R,
				aMatrix[componentIndex][componentIndex],
				bMatrix[componentIndex][componentIndex], concentration[i],
				T[i]);

	return CpkOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::hT(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::hkT(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::Fv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> FvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = RedlichKwongFluid::Fmx(hPlanck, molecMass, kB, T[i],
				rearX[i], molecVolumeMixture[i], aMatrixMolec, bMatrixMolec);

	return FvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::Fvk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> FvOutput(concentration.size());

	const auto & Mk = M[componentIndex];

	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < FvOutput.size(); ++i)
		FvOutput[i] = RedlichKwongFluid::Fv(concentration[i], T[i], R, Mk, ak,
				bk, hPlanck);

	return FvOutput;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::Sv(
		const std::vector<const std::valarray<scalar>*> & concentrations,
		const std::valarray<scalar> & T) const noexcept
{
	std::valarray<scalar> SvOutput(concentrations[0]->size());

	const auto X = calcMolarFrac(concentrations);

	const auto molecConc = (*concentrations[0]) * NAvogardro;

	const auto molecVolumeMixture = 1 / molecConc;

	const auto rearX = rearrangeMolFrac(X);

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = RedlichKwongFluid::Smx(hPlanck, molecMass, kB, T[i],
				rearX[i], molecVolumeMixture[i], aMatrixMolec, bMatrixMolec);

	return SvOutput * molecConc;
}

std::valarray<schemi::scalar> schemi::mixtureRedlichKwong::Svk(
		const std::valarray<scalar> & concentration,
		const std::valarray<scalar> & T,
		const std::size_t componentIndex) const noexcept
{
	std::valarray<scalar> SvOutput(concentration.size());

	const auto & Mk = M[componentIndex];

	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	for (std::size_t i = 0; i < SvOutput.size(); ++i)
		SvOutput[i] = RedlichKwongFluid::Sv(concentration[i], T[i], R, Mk, ak,
				bk, hPlanck);

	return SvOutput;
}

schemi::scalar schemi::mixtureRedlichKwong::Cv(
		const std::valarray<scalar> & concentrations) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	scalar CvMixOutput { 0 };

	for (std::size_t k = 0; k < X.size(); ++k)
		CvMixOutput += X[k] * CvArr[k];

	return CvMixOutput;
}

schemi::scalar schemi::mixtureRedlichKwong::pFromUv(
		const std::valarray<scalar> & concentrations, const scalar Uv) const
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

	return RedlichKwongFluid::pFromUv(R, CvMixture, Uv, concentrations[0],
			aMixture, bMixture);
}

schemi::scalar schemi::mixtureRedlichKwong::UvFromp(
		const std::valarray<scalar> & concentrations, const scalar p) const
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

	return RedlichKwongFluid::UvFromp(R, CvMixture, p, concentrations[0],
			aMixture, bMixture);
}

schemi::scalar schemi::mixtureRedlichKwong::pcFromT(
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

	return RedlichKwongFluid::pcFromT(R, T, concentrations[0], aMixture,
			bMixture);
}

schemi::scalar schemi::mixtureRedlichKwong::pcFromTk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	return RedlichKwongFluid::pcFromT(R, T, concentration, ak, bk);
}

schemi::scalar schemi::mixtureRedlichKwong::UvcFromT(
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

	return RedlichKwongFluid::UvcFromT(CvMixture, T, concentrations[0],
			aMixture, bMixture);
}

schemi::scalar schemi::mixtureRedlichKwong::UvcFromTk(
		const scalar concentration, const scalar T,
		const std::size_t componentIndex) const noexcept
{
	const auto & Cvk = CvArr[componentIndex];

	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	return RedlichKwongFluid::UvcFromT(Cvk, T, concentration, ak, bk);
}

schemi::scalar schemi::mixtureRedlichKwong::TFromUv(
		const std::valarray<scalar> & concentrations, const scalar Uv) const
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

	return RedlichKwongFluid::TFromUv(CvMixture, Uv, concentrations[0],
			aMixture, bMixture);
}

schemi::scalar schemi::mixtureRedlichKwong::dpdrho(
		const std::valarray<scalar> & concentrations, const scalar Uv) const
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

	return RedlichKwongFluid::dpdrho(CvMixture, R / CvMixture, Uv,
			concentrations[0], MMixture, aMixture, bMixture);
}

schemi::scalar schemi::mixtureRedlichKwong::dpdUv(
		const std::valarray<scalar> & concentrations, const scalar Uv) const
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

	return RedlichKwongFluid::dpdUv(CvMixture, R / CvMixture, concentrations[0],
			Uv, aMixture, bMixture);
}

schemi::scalar schemi::mixtureRedlichKwong::nonIdeality(
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

	return RedlichKwongFluid::nonIdeality(concentrations[0], aMixture, bMixture,
			T);
}

schemi::scalar schemi::mixtureRedlichKwong::sqSonicSpeed(
		const std::valarray<scalar> & concentrations, const scalar density,
		const scalar Uv, const scalar pressure) const noexcept
{
	return dpdrho(concentrations, Uv)
			+ dpdUv(concentrations, Uv) * (Uv + pressure) / density;
}

schemi::scalar schemi::mixtureRedlichKwong::Cp(
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

	return RedlichKwongFluid::Cp(CvMixture, R, aMixture, bMixture,
			concentrations[0], T);
}

schemi::scalar schemi::mixtureRedlichKwong::Cpk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	return RedlichKwongFluid::Cp(CvArr[componentIndex], R,
			aMatrix[componentIndex][componentIndex],
			bMatrix[componentIndex][componentIndex], concentration, T);
}

schemi::scalar schemi::mixtureRedlichKwong::hT(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	return (UvcFromT(concentrations, T) + pcFromT(concentrations, T)) / T;
}

schemi::scalar schemi::mixtureRedlichKwong::hkT(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	return (UvcFromTk(concentration, T, componentIndex)
			+ pcFromTk(concentration, T, componentIndex)) / T;
}

schemi::scalar schemi::mixtureRedlichKwong::Fv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	return RedlichKwongFluid::Fmx(hPlanck, molecMass, kB, T, X,
			molecVolumeMixture, aMatrixMolec, bMatrixMolec) * molecConc;
}

schemi::scalar schemi::mixtureRedlichKwong::Fvk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	return RedlichKwongFluid::Fv(concentration, T, R, M[componentIndex], ak, bk,
			hPlanck);
}

schemi::scalar schemi::mixtureRedlichKwong::Sv(
		const std::valarray<scalar> & concentrations,
		const scalar T) const noexcept
{
	const auto X = calcMolarFrac(concentrations);

	const auto molecConc { concentrations[0] * NAvogardro };

	const auto molecVolumeMixture { 1 / molecConc };

	return RedlichKwongFluid::Smx(hPlanck, molecMass, kB, T, X,
			molecVolumeMixture, aMatrixMolec, bMatrixMolec) * molecConc;
}

schemi::scalar schemi::mixtureRedlichKwong::Svk(const scalar concentration,
		const scalar T, const std::size_t componentIndex) const noexcept
{
	const auto & ak = aMatrix[componentIndex][componentIndex];
	const auto & bk = bMatrix[componentIndex][componentIndex];

	return RedlichKwongFluid::Sv(concentration, T, R, M[componentIndex], ak, bk,
			hPlanck);
}
