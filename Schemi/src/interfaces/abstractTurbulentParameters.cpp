/*
 * abstractTurbulentParameters.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractTurbulentParameters.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include "errorsEnum.hpp"
#include "intExpPow.hpp"
#include "exception.hpp"
#include "vectorVectorDotProduct.hpp"

schemi::scalar schemi::abstractTurbulentParameters::Cmu() const noexcept
{
	return Cmu_;
}
schemi::scalar schemi::abstractTurbulentParameters::C0() const noexcept
{
	return C0_;
}
schemi::scalar schemi::abstractTurbulentParameters::C1() const noexcept
{
	return C1_;
}
schemi::scalar schemi::abstractTurbulentParameters::C2() const noexcept
{
	return C2_;
}
schemi::scalar schemi::abstractTurbulentParameters::C3() const noexcept
{
	return C3_;
}
schemi::scalar schemi::abstractTurbulentParameters::sigmaSc() const noexcept
{
	return sigmaSc_;
}
schemi::scalar schemi::abstractTurbulentParameters::sigmaT() const noexcept
{
	return sigmaT_;
}
schemi::scalar schemi::abstractTurbulentParameters::sigmaE() const noexcept
{
	return sigmaE_;
}
schemi::scalar schemi::abstractTurbulentParameters::sigmaK() const noexcept
{
	return sigmaK_;
}
schemi::scalar schemi::abstractTurbulentParameters::sigmaEps() const noexcept
{
	return sigmaEps_;
}
schemi::scalar schemi::abstractTurbulentParameters::sigmaa() const noexcept
{
	return sigmaa_;
}
schemi::scalar schemi::abstractTurbulentParameters::sigmab() const noexcept
{
	return sigmab_;
}
schemi::scalar schemi::abstractTurbulentParameters::Ca1() const noexcept
{
	return Ca1_;
}
schemi::scalar schemi::abstractTurbulentParameters::Cb1() const noexcept
{
	return Cb1_;
}
schemi::scalar schemi::abstractTurbulentParameters::mink() const noexcept
{
	return mink_;
}
schemi::scalar schemi::abstractTurbulentParameters::mineps() const noexcept
{
	return mineps_;
}
schemi::scalar schemi::abstractTurbulentParameters::CMS_R() const noexcept
{
	return CMS_R_;
}
schemi::scalar schemi::abstractTurbulentParameters::CMS_D() const noexcept
{
	return CMS_D_;
}
schemi::scalar schemi::abstractTurbulentParameters::CMS_B() const noexcept
{
	return CMS_B_;
}

schemi::abstractTurbulentParameters::abstractTurbulentParameters(

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

const scalar CMS_B_In) noexcept :

		Cmu_ { CmuI },

		C0_ { C0In },

		C1_ { C1In },

		C2_ { C2In },

		C3_ { C3In },

		sigmaSc_ { sigmaScIn },

		sigmaT_ { sigmaTIn },

		sigmaE_ { sigmaEIn },

		sigmaK_ { sigmakIn },

		sigmaEps_ { sigmaepsIn },

		sigmaa_ { sigmaaIn },

		sigmab_ { sigmabIn },

		Ca1_ { Ca1In },

		Cb1_ { Cb1In },

		mink_ { minkIn },

		mineps_ { mienpsIn },

		CMS_R_ { CMS_R_In },

		CMS_D_ { CMS_D_In },

		CMS_B_ { CMS_B_In }
{
}

void schemi::abstractTurbulentParameters::readTurbulentParameters(
		std::ifstream & turbulentParametersFile)
{
	std::string skipBuffer;

	std::string thetaBType;

	turbulentParametersFile >> skipBuffer >> Cmu_ >> skipBuffer >> C0_
			>> skipBuffer >> C1_ >> skipBuffer >> C2_ >> skipBuffer >> C3_
			>> skipBuffer >> sigmaSc_ >> skipBuffer >> sigmaT_ >> skipBuffer
			>> sigmaE_ >> skipBuffer >> sigmaK_ >> skipBuffer >> sigmaEps_
			>> skipBuffer >> sigmaa_ >> skipBuffer >> sigmab_ >> skipBuffer
			>> Ca1_ >> skipBuffer >> Cb1_ >> skipBuffer >> mink_ >> skipBuffer
			>> mineps_ >> skipBuffer >> CMS_R_ >> skipBuffer >> CMS_D_
			>> skipBuffer >> CMS_B_ >> skipBuffer >> thetaBType;

	if (thetaBType == "sqrt")
	{
		thetaB_pointer = [this](const vector & a, const scalar k,
				const scalar /*epsilon*/, const vector& /*gradMav_n*/,
				const scalar /*a_s2*/,
				const std::pair<scalar, vector>& /*rho, gradRho*/,
				const std::pair<scalar, vector>& /*p, gradP*/,
				const scalar /*nu_t*/) 
				{
					return 1. / std::sqrt(1 + (a & a) / (CMS_B() * k));
				};

		std::cout << "thetaB is calculated with sqrt function." << std::endl;
	}
	else if (thetaBType == "cbrt")
	{
		thetaB_pointer = [this](const vector & a, const scalar k,
				const scalar /*epsilon*/, const vector& /*gradMav_n*/,
				const scalar /*a_s2*/,
				const std::pair<scalar, vector>& /*rho, gradRho*/,
				const std::pair<scalar, vector>& /*p, gradP*/,
				const scalar /*nu_t*/) 
				{
					return 1. / std::cbrt(1 + (a & a) / (CMS_B() * k));
				};

		std::cout << "thetaB is calculated with cbrt function." << std::endl;
	}
	else if (thetaBType == "gradMav")
	{
		thetaB_pointer =
				[this](const vector& /*a*/, const scalar k,
						const scalar epsilon, const vector & gradMav_n,
						const scalar /*a_s2*/,
						const std::pair<scalar, vector>& /*rho, gradRho*/,
						const std::pair<scalar, vector>& /*p, gradP*/,
						const scalar /*nu_t*/) 
						{
							const auto k3eps2 = pow<scalar, 3>(k) / pow<scalar, 2>(epsilon);
							return 1. / std::sqrt(1 + (gradMav_n & gradMav_n) * k3eps2 / CMS_B());
						};

		std::cout
				<< "thetaB is calculated with gradient of average molar mass function."
				<< std::endl;
	}
	else if (thetaBType == "wD")
	{
		thetaB_pointer =
				[this](const vector& /*a*/, const scalar k,
						const scalar /*epsilon*/, const vector& /*gradMav_n*/,
						const scalar a_s2,
						const std::pair<scalar, vector> & rhoGradRho,
						const std::pair<scalar, vector> & pGradP,
						const scalar nu_t) 
						{
							const auto wD = -nu_t * (rhoGradRho.second/rhoGradRho.first -pGradP.second/(rhoGradRho.first * a_s2));

							return 1. / std::sqrt(1 + (wD & wD) / (CMS_B()*k));
						};

		std::cout
				<< "thetaB is calculated with algebraically calculated turbulent mass flux."
				<< std::endl;
	}
	else if (thetaBType == "gradMavGradRho")
	{
		thetaB_pointer =
				[this](const vector& /*a*/, const scalar k,
						const scalar epsilon, const vector & gradMav_n,
						const scalar /*a_s2*/,
						const std::pair<scalar, vector> & rhoGradRho,
						const std::pair<scalar, vector> & pGradP,
						const scalar /*nu_t*/) 
						{
							const auto k4eps2 = pow<scalar, 4>(k) / pow<scalar, 2>(epsilon);
							return 1. / std::sqrt(1 + std::abs(gradMav_n & rhoGradRho.second)/pGradP.first * k4eps2 / CMS_B());
						};

		std::cout
				<< "thetaB is calculated with gradient of average molar mass and density gradient function."
				<< std::endl;
	}
	else if (thetaBType == "no")
		std::cout << "thetaB is 1." << std::endl;
	else
		throw exception(
				"Wrong type of thetaB function. Must be <<sqrt>> or <<cbrt>> or <<gradMav>> or <<wD>> or <<gradMavGradRho>> or <<no>>.",
				errorsEnum::initializationError);
}

schemi::scalar schemi::abstractTurbulentParameters::thetaS(
		const scalar divV) const noexcept
{
	return (divV > -CMS_R());
}

schemi::scalar schemi::abstractTurbulentParameters::thetaB(const vector & a,
		const scalar k, const scalar epsilon, const vector & gradMav_n,
		const scalar a_s2, const std::pair<scalar, vector> & rhoGradRho,
		const std::pair<scalar, vector> & pGradP,
		const scalar nu_t) const noexcept
{
	return thetaB_pointer(a, k, epsilon, gradMav_n, a_s2, rhoGradRho, pGradP,
			nu_t);
}

schemi::abstractTurbulentParameters::~abstractTurbulentParameters() noexcept
{
}
