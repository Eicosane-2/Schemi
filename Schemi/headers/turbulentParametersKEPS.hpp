/*
 * turbulentParametersKEPS.hpp
 *
 *  Created on: 2019/12/06
 *      Author: Maxim Boldyrev
 *
 *      Class for k-epsilon family models parameters.
 */

#ifndef TURBULENTPARAMETERSKEPS_HPP_
#define TURBULENTPARAMETERSKEPS_HPP_

#include "abstractTurbulentParameters.hpp"

namespace schemi
{
class turbulentParametersKEPS: public abstractTurbulentParameters
{
	scalar thetaS(const scalar divV, const scalar k, const scalar eps,
			const scalar CMS_par) const noexcept;
public:
	explicit turbulentParametersKEPS(

	const scalar CmuI = 0.09,

	const scalar C0In = 0.95,

	const scalar C1In = 1.44,

	const scalar C2In = 1.92,

	const scalar C3In = 2.0,

	const scalar sigmaScIn = 0.5,

	const scalar sigmaTIn = 1,

	const scalar sigmaEIn = 1,

	const scalar sigmakIn = 0.5,

	const scalar sigmaepsIn = 0.4,

	const scalar sigmaaIn = 0.2,

	const scalar sigmabIn = 0.5,

	const scalar Ca1In = 3.0,

	const scalar Cb1In = 3.0,

	const scalar minkIn = 1E-15,

	const scalar mienpsIn = 1E-15,

	const scalar CMS_R_In = 0.1,

	const scalar CMS_D_In = 1.0,

	const scalar CMS_B_In = 0.01) noexcept;

	std::valarray<scalar> calculateNut(const std::valarray<scalar> & k,
			const std::valarray<scalar> & eps) const noexcept override;

	scalar calculateNut(const scalar k, const scalar eps) const noexcept
			override;

	std::valarray<scalar> rhoepsilon(
			const bunchOfFields<cubicCell> & cf) const noexcept override;

	scalar thetaS_R(const scalar divV, const scalar k,
			const scalar eps) const noexcept override;

	volumeField<scalar> thetaS_R(const volumeField<scalar> & divV,
			const volumeField<scalar> & k,
			const volumeField<scalar> & eps) const noexcept override;

	surfaceField<scalar> thetaS_R(const surfaceField<scalar> & divV,
			const surfaceField<scalar> & k,
			const surfaceField<scalar> & eps) const noexcept override;

	scalar thetaS_D(const scalar divV, const scalar k,
			const scalar eps) const noexcept override;

	volumeField<scalar> thetaS_D(const volumeField<scalar> & divV,
			const volumeField<scalar> & k,
			const volumeField<scalar> & eps) const noexcept override;

	surfaceField<scalar> thetaS_D(const surfaceField<scalar> & divV,
			const surfaceField<scalar> & k,
			const surfaceField<scalar> & eps) const noexcept override;
};
}  // namespace schemi

#endif /* TURBULENTPARAMETERSKEPS_HPP_ */
