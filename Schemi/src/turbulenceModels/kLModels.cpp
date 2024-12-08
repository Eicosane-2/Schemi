/*
 * kLModels.cpp
 *
 *  Created on: 2024/12/03
 *      Author: Maxim Boldyrev
 */

#include "kLModels.hpp"

#include "volumeField.hpp"
#include "surfaceField.hpp"
#include "intExpPow.hpp"

schemi::scalar schemi::kLModels::thetaS(const scalar divV, const scalar k,
		const scalar L, const scalar CMS) const noexcept
{
	const auto Cmu2 = Cmu() * Cmu();
	const auto k2 = k * k;
	const auto eps2 = pow<scalar, 3>(k) / pow<scalar, 2>(L);

	const auto CMS2 = CMS * CMS;

	const auto divV2 = divV * divV;

	return 1. / std::sqrt(1 + 16. / 9. * Cmu2 * k2 / (CMS2 * eps2) * divV2);
}

void schemi::kLModels::CMS_RSet(const scalar v) noexcept
{
	CMS_R_ = v;
}

void schemi::kLModels::CMS_DSet(const scalar v) noexcept
{
	CMS_D_ = v;
}

schemi::scalar schemi::kLModels::CMS_R() const noexcept
{
	return CMS_R_;
}

schemi::scalar schemi::kLModels::CMS_D() const noexcept
{
	return CMS_D_;
}

schemi::kLModels::kLModels(const mesh & meshIn, const bool turb, const bool a,
		const bool b, const turbulenceModel model, const scalar Cmu_in,
		const scalar sigmaSc_in, const scalar sigmaT_in, const scalar sigmaE_in,
		const scalar sigmak_in, const scalar sigmaL_in, const scalar sigmaa_in,
		const scalar sigmab_in, const scalar CMS_R_in,
		const scalar CMS_D_in) noexcept :
		abstractTurbulenceModel(meshIn, turb, a, b, model, Cmu_in, sigmaSc_in,
				sigmaT_in, sigmaE_in, sigmak_in, sigmaL_in, sigmaa_in,
				sigmab_in), CMS_R_(CMS_R_in), CMS_D_(CMS_D_in)
{
}

std::valarray<schemi::scalar> schemi::kLModels::calculateNut(
		const std::valarray<scalar> & k,
		const std::valarray<scalar> & L) const noexcept
{
	return Cmu() * std::sqrt(k) * L;
}

schemi::scalar schemi::kLModels::calculateNut(const scalar k,
		const scalar L) const noexcept
{
	return Cmu() * std::sqrt(k) * L;
}

std::valarray<schemi::scalar> schemi::kLModels::rhoepsilon(
		const bunchOfFields<cubicCell> & cf) const noexcept
{
	return cf.density[0]() * std::sqrt(cf.kTurb() * cf.kTurb() * cf.kTurb())
			/ cf.epsTurb();
}

schemi::scalar schemi::kLModels::thetaS_R(const scalar divV, const scalar k,
		const scalar L) const noexcept
{
	return thetaS(divV, k, L, CMS_R());
}

schemi::volumeField<schemi::scalar> schemi::kLModels::thetaS_R(
		const volumeField<scalar> & divV, const volumeField<scalar> & k,
		const volumeField<scalar> & L) const noexcept
{
	volumeField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		returnField.r()[i] = thetaS(divV()[i], k()[i], L()[i], CMS_R());

	return returnField;
}

schemi::surfaceField<schemi::scalar> schemi::kLModels::thetaS_R(
		const surfaceField<scalar> & divV, const surfaceField<scalar> & k,
		const surfaceField<scalar> & L) const noexcept
{
	surfaceField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		returnField.r()[i] = thetaS(divV()[i], k()[i], L()[i], CMS_R());

	return returnField;
}

schemi::scalar schemi::kLModels::thetaS_D(const scalar divV, const scalar k,
		const scalar L) const noexcept
{
	if (divV < 0.)
		return thetaS(divV, k, L, CMS_D());
	else
		return 1.;
}

schemi::volumeField<schemi::scalar> schemi::kLModels::thetaS_D(
		const volumeField<scalar> & divV, const volumeField<scalar> & k,
		const volumeField<scalar> & L) const noexcept
{
	volumeField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		if (divV()[i] < 0.)
			returnField.r()[i] = thetaS(divV()[i], k()[i], L()[i], CMS_D());
		else
			returnField.r()[i] = 1.;

	return returnField;
}

schemi::surfaceField<schemi::scalar> schemi::kLModels::thetaS_D(
		const surfaceField<scalar> & divV, const surfaceField<scalar> & k,
		const surfaceField<scalar> & L) const noexcept
{
	surfaceField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		if (divV()[i] < 0.)
			returnField.r()[i] = thetaS(divV()[i], k()[i], L()[i], CMS_D());
		else
			returnField.r()[i] = 1.;

	return returnField;
}

schemi::kLModels::~kLModels() noexcept
{
}
