/*
 * turbulentParametersKEPS.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "turbulentParametersKEPS.hpp"

schemi::scalar schemi::turbulentParametersKEPS::thetaS(const scalar divV,
		const scalar k, const scalar eps, const scalar CMS_par) const noexcept
{
	const auto Cmu2 = Cmu() * Cmu();
	const auto k2 = k * k;
	const auto eps2 = eps * eps;

	const auto CMS2 = CMS_par * CMS_par;

	const auto divV2 = divV * divV;

	return 1. / std::sqrt(1 + 16. / 9. * Cmu2 * k2 / (CMS2 * eps2) * divV2);
}

schemi::turbulentParametersKEPS::turbulentParametersKEPS(

const mesh & meshIn,

const scalar CmuIn,

const scalar C0In,

const scalar C1In,

const scalar C2In,

const scalar C3In,

const scalar sigmaScIn,

const scalar sigmaTIn,

const scalar sigmaEIn,

const scalar sigmakIn,

const scalar sigmaepsIn,

const scalar sigmaaIn,

const scalar sigmabIn,

const scalar Ca1In,

const scalar Cb1In,

const scalar minkIn,

const scalar mienpsIn,

const scalar CMS_R_In,

const scalar CMS_D_In,

const scalar CMS_B_In) noexcept :
		abstractTurbulentParameters(meshIn, CmuIn, C0In, C1In, C2In, C3In,
				sigmaScIn, sigmaTIn, sigmaEIn, sigmakIn, sigmaepsIn, sigmaaIn,
				sigmabIn, Ca1In, Cb1In, minkIn, mienpsIn, CMS_R_In, CMS_D_In,
				CMS_B_In)
{
}

std::valarray<schemi::scalar> schemi::turbulentParametersKEPS::calculateNut(
		const std::valarray<scalar> & k,
		const std::valarray<scalar> & eps) const noexcept
{
	return Cmu() * k * k / eps;
}

schemi::scalar schemi::turbulentParametersKEPS::calculateNut(const scalar k,
		const scalar eps) const noexcept
{
	return Cmu() * k * k / eps;
}

std::valarray<schemi::scalar> schemi::turbulentParametersKEPS::rhoepsilon(
		const bunchOfFields<cubicCell> & cf) const noexcept
{
	return cf.rhoepsTurb();
}

schemi::scalar schemi::turbulentParametersKEPS::thetaS_R(const scalar divV,
		const scalar k, const scalar eps) const noexcept
{
	return thetaS(divV, k, eps, CMS_R());
}

schemi::volumeField<schemi::scalar> schemi::turbulentParametersKEPS::thetaS_R(
		const volumeField<scalar> & divV, const volumeField<scalar> & k,
		const volumeField<scalar> & eps) const noexcept
{
	volumeField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		returnField.r()[i] = thetaS(divV()[i], k()[i], eps()[i], CMS_R());

	return returnField;
}

schemi::surfaceField<schemi::scalar> schemi::turbulentParametersKEPS::thetaS_R(
		const surfaceField<scalar> & divV, const surfaceField<scalar> & k,
		const surfaceField<scalar> & eps) const noexcept
{
	surfaceField<scalar> returnField { divV.meshRef(), 0. };

	for (std::size_t i = 0; i < returnField.size(); ++i)
		returnField.r()[i] = thetaS(divV()[i], k()[i], eps()[i], CMS_R());

	return returnField;
}

schemi::scalar schemi::turbulentParametersKEPS::thetaS_D(const scalar divV,
		const scalar k, const scalar eps) const noexcept
{
	if (divV < 0.)
		return thetaS(divV, k, eps, CMS_D());
	else
		return 1.;
}

schemi::volumeField<schemi::scalar> schemi::turbulentParametersKEPS::thetaS_D(
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

schemi::surfaceField<schemi::scalar> schemi::turbulentParametersKEPS::thetaS_D(
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
