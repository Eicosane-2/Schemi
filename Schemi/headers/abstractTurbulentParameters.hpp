/*
 * abstractTurbulentParameters.hpp
 *
 *  Created on: 2021/05/06
 *      Author: Maxim Boldyrev
 *
 *      Interface class for turbulence model parameters and fields.
 */

#ifndef ABSTRACTTURBULENTPARAMETERS_HPP_
#define ABSTRACTTURBULENTPARAMETERS_HPP_

#include <functional>
#include <utility>
#include <fstream>

#include "cubicCell.hpp"
#include "bunchOfFields.hpp"
#include "surfaceField.hpp"
#include "volumeField.hpp"

namespace schemi
{
class boundaryConditionValue;

class abstractTurbulentParameters
{
	std::function<
			scalar(const vector& /*a*/, const scalar /*k*/,
					const scalar /*epsilon*/, const vector& /*gradMav_n*/,
					const scalar /*a_s2*/,
					const std::pair<scalar, vector>&& /*rho, gradRho*/,
					const std::pair<scalar, vector>&& /*p, gradP*/,
					const scalar /*nu_t*/)> thetaB_pointer = [this](
			const vector&, const scalar, const scalar, const vector&,
			const scalar, const std::pair<scalar, vector>&&,
			const std::pair<scalar, vector>&&, const scalar) 
			{
				return 1.;
			};

	scalar Cmu_;
	scalar C0_;
	scalar C1_;
	scalar C2_;
	scalar C3_;
	scalar sigmaSc_;
	scalar sigmaT_;
	scalar sigmaE_;
	scalar sigmaK_;
	scalar sigmaEps_;
	scalar sigmaa_;
	scalar sigmab_;
	scalar Ca1_;
	scalar Cb1_;
	scalar mink_;
	scalar mineps_;
	scalar CMS_R_;
	scalar CMS_D_;
	scalar CMS_B_;

	volumeField<scalar> y;

	std::valarray<scalar> solveMatrix(const volumeField<scalar> & phi,
			const boundaryConditionValue & boundCond) const;

public:
	scalar Cmu() const noexcept;
	scalar C0() const noexcept;
	scalar C1() const noexcept;
	scalar C2() const noexcept;
	scalar C3() const noexcept;
	scalar sigmaSc() const noexcept;
	scalar sigmaT() const noexcept;
	scalar sigmaE() const noexcept;
	scalar sigmaK() const noexcept;
	scalar sigmaEps() const noexcept;
	scalar sigmaa() const noexcept;
	scalar sigmab() const noexcept;
	scalar Ca1() const noexcept;
	scalar Cb1() const noexcept;
	scalar mink() const noexcept;
	scalar mineps() const noexcept;
	scalar CMS_R() const noexcept;
	scalar CMS_D() const noexcept;
	scalar CMS_B() const noexcept;

	virtual ~abstractTurbulentParameters() noexcept =0;

	abstractTurbulentParameters(

	const mesh & meshIn,

	const scalar CmuI,

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

	const scalar CMS_B_In) noexcept;

	virtual std::valarray<scalar> calculateNut(
			const std::valarray<scalar>& /*k*/,
			const std::valarray<scalar>& /*eps*/) const noexcept =0;

	virtual scalar calculateNut(const scalar /*k*/,
			const scalar /*eps*/) const noexcept =0;

	void readTurbulentParameters(std::ifstream & turbulentParametersFile);

	virtual std::valarray<scalar> rhoepsilon(
			const bunchOfFields<cubicCell>& /*cf*/) const noexcept =0;

	virtual scalar thetaS_R(const scalar /*divV*/, const scalar /*k*/,
			const scalar /*eps*/) const noexcept =0;

	virtual volumeField<scalar> thetaS_R(const volumeField<scalar>& /*divV*/,
			const volumeField<scalar>& /*k*/,
			const volumeField<scalar>& /*eps*/) const noexcept =0;

	virtual surfaceField<scalar> thetaS_R(const surfaceField<scalar>& /*divV*/,
			const surfaceField<scalar>& /*k*/,
			const surfaceField<scalar>& /*eps*/) const noexcept =0;

	virtual scalar thetaS_D(const scalar /*divV*/, const scalar /*k*/,
			const scalar /*eps*/) const noexcept =0;

	virtual volumeField<scalar> thetaS_D(const volumeField<scalar>& /*divV*/,
			const volumeField<scalar>& /*k*/,
			const volumeField<scalar>& /*eps*/) const noexcept =0;

	virtual surfaceField<scalar> thetaS_D(const surfaceField<scalar>& /*divV*/,
			const surfaceField<scalar>& /*k*/,
			const surfaceField<scalar>& /*eps*/) const noexcept =0;

	scalar thetaS(const scalar divV) const noexcept;

	template<typename typeOfField>
	field<typeOfField, scalar> thetaS(
			const field<typeOfField, scalar> & divV) const noexcept
	{
		field<typeOfField, scalar> returnField { divV.meshP(), 0. };

		std::transform(std::begin(divV()), std::end(divV()),
				std::begin(returnField.r()), [this](const auto div) 
				{	return this->thetaS(div);});

		return returnField;
	}

	scalar thetaB(const vector & a, const scalar k, const scalar epsilon,
			const vector & gradMav_n, const scalar a_s2,
			const std::pair<scalar, vector> && rhoGradRho,
			const std::pair<scalar, vector> && pGradP,
			const scalar nu_t) const noexcept;

	template<typename typeOfField>
	field<scalar, typeOfField> thetaB(const field<vector, typeOfField> & a,
			const field<scalar, typeOfField> & k,
			const field<scalar, typeOfField> & epsilon,
			const field<vector, typeOfField> & gradMav_n,
			const std::valarray<scalar> & a_s2,
			const std::pair<field<scalar, typeOfField>,
					field<vector, typeOfField>> & rhoGradRho,
			const std::pair<field<scalar, typeOfField>,
					field<vector, typeOfField>> & pGradP,
			const field<scalar, typeOfField> & nu_t) const noexcept
	{
		field<scalar, typeOfField> returnField { a.meshP(), 0. };

		for (std::size_t i = 0; i < returnField.size(); ++i)
			returnField.r()[i] = thetaB(a()[i], k()[i], epsilon()[i],
					gradMav_n()[i], a_s2[i],
					std::pair<scalar, vector>(rhoGradRho.first()[i],
							rhoGradRho.second()[i]),
					std::pair<scalar, vector>(pGradP.first()[i],
							pGradP.second()[i]), nu_t()[i]);

		return returnField;
	}

	constexpr static scalar minb_value { 1E-15 };

	void calculateNearWallDistance(const volumeField<scalar> & eps,
			const boundaryConditionValue & boundCond);
};
}  // namespace schemi

#endif /* ABSTRACTTURBULENTPARAMETERS_HPP_ */
