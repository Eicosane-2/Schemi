/*
 * turbulentParametersKL.hpp
 *
 *  Created on: 2019/12/06
 *      Author: Maxim Boldyrev
 *
 *      Class for k-L family models parameters.
 */

#ifndef TURBULENTPARAMETERSKL_HPP_
#define TURBULENTPARAMETERSKL_HPP_

#include "abstractTurbulentParameters.hpp"

namespace schemi
{
class turbulentParametersKL: public abstractTurbulentParameters
{
	scalar thetaS(const scalar divV, const scalar k, const scalar eps,
			const scalar CMS_par) const noexcept;
public:
	constexpr static scalar Dr = 3.11;

	explicit turbulentParametersKL(

	const mesh & meshIn,

	const scalar CmuI = 0.28,

	const scalar C0In = 1.5 - 1.05,

	const scalar C1In = 1.5 - 1.2,

	const scalar C2In = 1.5 - 1.92,

	const scalar C3In = 1.5 - 1.67,

	const scalar sigmaScIn = 0.5,

	const scalar sigmaTIn = 1,

	const scalar sigmaEIn = 1,

	const scalar sigmakIn = 1.0,

	const scalar sigmaepsIn = 0.1,

	const scalar sigmaaIn = 1.0,

	const scalar sigmabIn = 1.0,

	const scalar Ca1In = 1.9,

	const scalar Ca1maxIn = 19.0,

	const scalar Cb1In = 1.2,

	const scalar minkIn = 1E-15,

	const scalar mienpsIn = 1E-15,

	const scalar CMS_R_In = 0.1,

	const scalar CMS_D_In = 1.0,

	const scalar CMS_A_In = 1.0,

	const scalar CMS_M_In = 1.0) noexcept;

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

#endif /* TURBULENTPARAMETERSKL_HPP_ */
