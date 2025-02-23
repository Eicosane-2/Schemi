/*
 * kEpsModels.cpp
 *
 *  Created on: 2024/12/03
 *      Author: Maxim Boldyrev
 */

#include "kEpsModels.hpp"

#include "volumeField.hpp"
#include "surfaceField.hpp"

schemi::scalar schemi::kEpsModels::thetaS(const scalar divV, const scalar k,
		const scalar eps, const scalar CMS) const noexcept
{
	const auto Cmu2 = Cmu() * Cmu();
	const auto k2 = k * k;
	const auto eps2 = eps * eps;

	const auto CMS2 = CMS * CMS;

	const auto divV2 = divV * divV;

	return 1. / std::sqrt(1 + 16. / 9. * Cmu2 * k2 / (CMS2 * eps2) * divV2);
}

void schemi::kEpsModels::CMS_RSet(const scalar v) noexcept
{
	CMS_R_ = v;
}

void schemi::kEpsModels::CMS_DSet(const scalar v) noexcept
{
	CMS_D_ = v;
}

schemi::scalar schemi::kEpsModels::CMS_R() const noexcept
{
	return CMS_R_;
}

schemi::scalar schemi::kEpsModels::CMS_D() const noexcept
{
	return CMS_D_;
}

schemi::kEpsModels::kEpsModels(const mesh & meshIn, const bool turb,
		const bool a, const bool b, const turbulenceModel model,
		const scalar Cmu_in, const scalar sigmaSc_in, const scalar sigmaT_in,
		const scalar sigmaE_in, const scalar sigmak_in,
		const scalar sigmaeps_in, const scalar sigmaa_in,
		const scalar sigmab_in, const scalar CMS_R_in,
		const scalar CMS_D_in) noexcept :
		abstractTurbulenceModel(meshIn, turb, a, b, model, Cmu_in, sigmaSc_in,
				sigmaT_in, sigmaE_in, sigmak_in, sigmaeps_in, sigmaa_in,
				sigmab_in), CMS_R_(CMS_R_in), CMS_D_(CMS_D_in)
{
}

std::valarray<schemi::scalar> schemi::kEpsModels::calculateNut(
		const std::valarray<scalar> & k,
		const std::valarray<scalar> & eps) const noexcept
{
	return Cmu() * k * k / eps;
}

schemi::scalar schemi::kEpsModels::calculateNut(const scalar k,
		const scalar eps) const noexcept
{
	return Cmu() * k * k / eps;
}

std::valarray<schemi::scalar> schemi::kEpsModels::rhoepsilon(
		const bunchOfFields<cubicCell> & cf,
		const abstractMixtureThermodynamics&) const noexcept
{
	return cf.rhoepsTurb();
}

schemi::scalar schemi::kEpsModels::thetaS_R(const scalar divV, const scalar k,
		const scalar eps) const noexcept
{
	return thetaS(divV, k, eps, CMS_R());
}

schemi::volumeField<schemi::scalar> schemi::kEpsModels::thetaS_R(
		const volumeField<scalar> & divV, const volumeField<scalar> & k,
		const volumeField<scalar> & eps) const noexcept
{
	volumeField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		returnField.r()[i] = thetaS(divV()[i], k()[i], eps()[i], CMS_R());

	return returnField;
}

schemi::surfaceField<schemi::scalar> schemi::kEpsModels::thetaS_R(
		const surfaceField<scalar> & divV, const surfaceField<scalar> & k,
		const surfaceField<scalar> & eps) const noexcept
{
	surfaceField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		returnField.r()[i] = thetaS(divV()[i], k()[i], eps()[i], CMS_R());

	return returnField;
}

schemi::scalar schemi::kEpsModels::thetaS_D(const scalar divV, const scalar k,
		const scalar eps) const noexcept
{
	if (divV < 0.)
		return thetaS(divV, k, eps, CMS_D());
	else
		return 1.;
}

schemi::volumeField<schemi::scalar> schemi::kEpsModels::thetaS_D(
		const volumeField<scalar> & divV, const volumeField<scalar> & k,
		const volumeField<scalar> & eps) const noexcept
{
	volumeField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		if (divV()[i] < 0.)
			returnField.r()[i] = thetaS(divV()[i], k()[i], eps()[i], CMS_D());
		else
			returnField.r()[i] = 1.;

	return returnField;
}

schemi::surfaceField<schemi::scalar> schemi::kEpsModels::thetaS_D(
		const surfaceField<scalar> & divV, const surfaceField<scalar> & k,
		const surfaceField<scalar> & eps) const noexcept
{
	surfaceField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		if (divV()[i] < 0.)
			returnField.r()[i] = thetaS(divV()[i], k()[i], eps()[i], CMS_D());
		else
			returnField.r()[i] = 1.;

	return returnField;
}

schemi::kEpsModels::~kEpsModels() noexcept
{
}
