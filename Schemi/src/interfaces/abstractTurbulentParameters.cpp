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

#include "boundaryConditionValue.hpp"
#include "errorsEnum.hpp"
#include "exception.hpp"
#include "vectorVectorDotProduct.hpp"
#include "intExpPow.hpp"
#include "gradient.hpp"
#include "SLEMatrix.hpp"
#include "conjugateGradientSolver.hpp"

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
schemi::scalar schemi::abstractTurbulentParameters::Ca1max() const noexcept
{
	return Ca1_max;
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
schemi::scalar schemi::abstractTurbulentParameters::CMS_A() const noexcept
{
	return CMS_A_;
}
schemi::scalar schemi::abstractTurbulentParameters::CMS_M() const noexcept
{
	return CMS_M_;
}

schemi::abstractTurbulentParameters::abstractTurbulentParameters(

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

const scalar Ca1maxIn,

const scalar Cb1In,

const scalar minkIn,

const scalar mienpsIn,

const scalar CMS_R_In,

const scalar CMS_D_In,

const scalar CMS_A_In,

const scalar CMS_M_In) noexcept :

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

		Ca1_max { Ca1maxIn },

		Cb1_ { Cb1In },

		mink_ { minkIn },

		mineps_ { mienpsIn },

		CMS_R_ { CMS_R_In },

		CMS_D_ { CMS_D_In },

		CMS_A_ { CMS_A_In },

		CMS_M_ { CMS_M_In },

		y { meshIn }
{
}

void schemi::abstractTurbulentParameters::readTurbulentParameters(
		std::ifstream & turbulentParametersFile)
{
	std::string skipBuffer;

	turbulentParametersFile >> skipBuffer >> Cmu_ >> skipBuffer >> C0_
			>> skipBuffer >> C1_ >> skipBuffer >> C2_ >> skipBuffer >> C3_
			>> skipBuffer >> sigmaSc_ >> skipBuffer >> sigmaT_ >> skipBuffer
			>> sigmaE_ >> skipBuffer >> sigmaK_ >> skipBuffer >> sigmaEps_
			>> skipBuffer >> sigmaa_ >> skipBuffer >> sigmab_ >> skipBuffer
			>> Ca1_ >> skipBuffer >> Ca1_max >> skipBuffer >> Cb1_ >> skipBuffer
			>> mink_ >> skipBuffer >> mineps_ >> skipBuffer >> CMS_R_
			>> skipBuffer >> CMS_D_ >> skipBuffer >> CMS_A_ >> skipBuffer
			>> CMS_M_;
}

schemi::scalar schemi::abstractTurbulentParameters::thetaS(
		const scalar divV) const noexcept
{
	return (divV > -CMS_R());
}

schemi::scalar schemi::abstractTurbulentParameters::thetaA(const vector & a,
		const scalar k, const scalar b) const noexcept
{
	const auto aa = (a & a);

	return aa / (2 * k) < CMS_A() && b < CMS_A() && aa / (2 * k) < b ?
			1 : Ca1max() / std::min(Ca1(), Cb1());
}

void schemi::abstractTurbulentParameters::calculateNearWallDistance(
		const volumeField<scalar> & eps,
		const boundaryConditionValue & boundCond)
{
	constexpr static scalar phiInit = 1;
	constexpr static scalar convergenceTolerance { convergenceToleranceGlobal };

	volumeField<scalar> phiOld = eps;

	phiOld.r() = phiInit;
	for (auto & b_i : phiOld.boundCond_r())
		if (b_i.first == boundaryConditionType::slip)
		{
			b_i.first = boundaryConditionType::fixedValueSurface;
			b_i.second = 0.;
		}
		else if ((b_i.first == boundaryConditionType::fixedValueCell)
				|| (b_i.first == boundaryConditionType::fixedValueSurface))
			b_i.first = boundaryConditionType::freeBoundary;

	boundCond.parallelism.correctBoundaryValues(phiOld);

	volumeField<scalar> phiNew = phiOld;

	std::size_t it { 0 };
	while (true)
	{
		it++;

		const scalar delta = std::abs((phiNew() - phiOld()) / phiOld()).max();

		if (((delta < convergenceTolerance) && (it > 1)) || (it > 10000))
			break;

		phiOld = phiNew;

		phiNew.r() = solveMatrix(phiOld, boundCond);

		boundCond.parallelism.correctBoundaryValues(phiNew);
	}

	const volumeField<vector> gradPhi = grad(phiNew, boundCond);

	volumeField<scalar> distance(phiNew.meshRef());

	for (std::size_t i = 0; i < distance.size(); ++i)
		distance.r()[i] = -gradPhi()[i].mag()
				+ std::sqrt(
						pow<scalar, 2>(gradPhi()[i].mag()) + 2 * phiNew()[i]);

	y = distance;
}

std::valarray<schemi::scalar> schemi::abstractTurbulentParameters::solveMatrix(
		const volumeField<scalar> & phi,
		const boundaryConditionValue & boundCond) const
{
	const auto & mesh_ = phi.meshRef();

	SLEMatrix wallDistanceMatrix(std::string("wallDistance"));

	wallDistanceMatrix.generateLaplacianSurfaceBoundary(phi,
			surfaceField<scalar>(mesh_, 1.0), boundCond);

	wallDistanceMatrix.SLE[0].freeTerm += 1;

	conjugateGradientSovler solv(100000, matrixSolver::conjugateGradient);

	return solv.solve(phi(), wallDistanceMatrix);
}

schemi::abstractTurbulentParameters::~abstractTurbulentParameters() noexcept
{
}
