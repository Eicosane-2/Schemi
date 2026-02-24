/*
 * abstractTurbulenceModel.cpp
 *
 *  Created on: 2024/12/03
 *      Author: Maxim Boldyrev
 */

#include "abstractTurbulenceModel.hpp"

#include <cmath>
#include <stdexcept>

#include "arithmeticAModel.hpp"
#include "BHRModel.hpp"
#include "BHRKLModel.hpp"
#include "decayModel.hpp"
#include "boundaryConditionValue.hpp"
#include "gradient.hpp"
#include "intExpPow.hpp"
#include "kEpsAModel.hpp"
#include "shearModel.hpp"
#include "SLEMatrix.hpp"
#include "conjugateGradientSolver.hpp"
#include "zeroModel.hpp"

void schemi::abstractTurbulenceModel::turbulenceSet(const bool v) noexcept
{
	turbulence_ = v;
}
void schemi::abstractTurbulenceModel::CmuSet(const scalar v) noexcept
{
	Cmu_ = v;
}
void schemi::abstractTurbulenceModel::sigmaScSet(const scalar v) noexcept
{
	sigmaSc_ = v;
}
void schemi::abstractTurbulenceModel::sigmaTSet(const scalar v) noexcept
{
	sigmaT_ = v;
}
void schemi::abstractTurbulenceModel::sigmaESet(const scalar v) noexcept
{
	sigmaE_ = v;
}
void schemi::abstractTurbulenceModel::sigmaKSet(const scalar v) noexcept
{
	sigmaK_ = v;
}
void schemi::abstractTurbulenceModel::sigmaEpsSet(const scalar v) noexcept
{
	sigmaEps_ = v;
}
void schemi::abstractTurbulenceModel::sigmaaSet(const scalar v) noexcept
{
	sigmaa_ = v;
}
void schemi::abstractTurbulenceModel::sigmabSet(const scalar v) noexcept
{
	sigmab_ = v;
}

std::valarray<schemi::scalar>& schemi::abstractTurbulenceModel::modelParametersSet() noexcept
{
	return modelParameters_;
}

const std::valarray<schemi::scalar>& schemi::abstractTurbulenceModel::modelParameters() const noexcept
{
	return modelParameters_;
}

void schemi::abstractTurbulenceModel::minkSet(const scalar v) noexcept
{
	mink_value_ = v;
}
void schemi::abstractTurbulenceModel::minepsilonSet(const scalar v) noexcept
{
	minepsilon_value_ = v;
}
void schemi::abstractTurbulenceModel::minbSet(const scalar v) noexcept
{
	minb_value_ = v;
}

bool schemi::abstractTurbulenceModel::turbulence() const noexcept
{
	return turbulence_;
}
bool schemi::abstractTurbulenceModel::aField() const noexcept
{
	return aField_;
}
bool schemi::abstractTurbulenceModel::bField() const noexcept
{
	return bField_;
}
schemi::turbulenceModel schemi::abstractTurbulenceModel::model() const noexcept
{
	return modelType_;
}
schemi::scalar schemi::abstractTurbulenceModel::Cmu() const noexcept
{
	return Cmu_;
}
schemi::scalar schemi::abstractTurbulenceModel::sigmaSc() const noexcept
{
	return sigmaSc_;
}
schemi::scalar schemi::abstractTurbulenceModel::sigmaT() const noexcept
{
	return sigmaT_;
}
schemi::scalar schemi::abstractTurbulenceModel::sigmaE() const noexcept
{
	return sigmaE_;
}
schemi::scalar schemi::abstractTurbulenceModel::sigmaK() const noexcept
{
	return sigmaK_;
}
schemi::scalar schemi::abstractTurbulenceModel::sigmaEps() const noexcept
{
	return sigmaEps_;
}
schemi::scalar schemi::abstractTurbulenceModel::sigmaa() const noexcept
{
	return sigmaa_;
}
schemi::scalar schemi::abstractTurbulenceModel::sigmab() const noexcept
{
	return sigmab_;
}
schemi::scalar schemi::abstractTurbulenceModel::mink() const noexcept
{
	return mink_value_;
}
schemi::scalar schemi::abstractTurbulenceModel::minepsilon() const noexcept
{
	return minepsilon_value_;
}
schemi::scalar schemi::abstractTurbulenceModel::minb() const noexcept
{
	return minb_value_;
}

schemi::abstractTurbulenceModel::abstractTurbulenceModel(const mesh & meshIn,
		const bool turb, const bool a, const bool b,
		const turbulenceModel model, const scalar Cmu_in,
		const scalar sigmaSc_in, const scalar sigmaT_in, const scalar sigmaE_in,
		const scalar sigmaK_in, const scalar sigmaEps_in,
		const scalar sigmaa_in, const scalar sigmab_in) noexcept :

		turbulence_(turb),

		aField_(a),

		bField_(b),

		modelType_(model),

		Cmu_ { Cmu_in },

		sigmaSc_ { sigmaSc_in },

		sigmaT_ { sigmaT_in },

		sigmaE_ { sigmaE_in },

		sigmaK_ { sigmaK_in },

		sigmaEps_ { sigmaEps_in },

		sigmaa_ { sigmaa_in },

		sigmab_ { sigmab_in },

		modelParameters_(0),

		y { meshIn }
{
}

void schemi::abstractTurbulenceModel::calculateNearWallDistance(
		const volumeField<scalar> & eps,
		const boundaryConditionValue & boundCond)
{
	constexpr static scalar phiInit = 1;
	constexpr static scalar convergenceTolerance { convergenceToleranceGlobal };

	volumeField<scalar> phiOld = eps;

	phiOld.val() = phiInit;
	for (auto & b_i : phiOld.boundCond_wr())
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

		const scalar delta = std::abs(
				(phiNew.cval() - phiOld.cval()) / phiOld.cval()).max();

		if (((delta < convergenceTolerance) && (it > 1)) || (it > 10000))
			break;

		phiOld = phiNew;

		phiNew.val() = solveMatrix(phiOld, boundCond);

		boundCond.parallelism.correctBoundaryValues(phiNew);
	}

	const volumeField<vector> gradPhi = grad(phiNew, boundCond);

	volumeField<scalar> distance(phiNew.meshRef());

	for (std::size_t i = 0; i < distance.size(); ++i)
		distance.val()[i] = -gradPhi.cval()[i].mag()
				+ std::sqrt(
						pow<scalar, 2>(gradPhi.cval()[i].mag())
								+ 2 * phiNew.cval()[i]);

	y = distance;
}

std::valarray<schemi::scalar> schemi::abstractTurbulenceModel::solveMatrix(
		const volumeField<scalar> & phi,
		const boundaryConditionValue & boundCond) const
{
	auto & mesh_ = phi.meshRef();

	SLEMatrix wallDistanceMatrix(std::string("wallDistance"));

	wallDistanceMatrix.generateLaplacianSurfaceBoundary(phi,
			surfaceField<scalar>(mesh_, 1.0), boundCond);

	wallDistanceMatrix.SLE[0].freeTerm += 1;

	conjugateGradientSovler solv(100000, matrixSolver::conjugateGradient);

	return solv.solve(phi.cval(), wallDistanceMatrix);
}

std::unique_ptr<schemi::abstractTurbulenceModel> schemi::abstractTurbulenceModel::createTurbulenceModel(
		const mesh & meshIn, const std::string & turbulenceONString,
		const std::string & sourceTypeString, const MPIHandler & parIn,
		const volumeField<vector> & uCellIn,
		const surfaceField<vector> & uSurfIn,
		const std::pair<std::size_t, std::string> & readDataPoint)
{
	bool turbulenceFlag;
	try
	{
		turbulenceFlag = onOffMap.at(turbulenceONString);
	} catch (std::out_of_range&)
	{
		throw exception("Unknown turbulence flag.",
				errors::initialisationError);
	}

	std::map<std::string, turbulenceModel> turbulenceModels;
	turbulenceModels.insert( { "zero", turbulenceModel::zeroModel });
	turbulenceModels.insert( { "decay", turbulenceModel::decayModel });
	turbulenceModels.insert( { "shear", turbulenceModel::shearModel });
	turbulenceModels.insert(
			{ "arithmetic1", turbulenceModel::arithmeticA1Model });
	turbulenceModels.insert(
			{ "arithmetic2", turbulenceModel::arithmeticA2Model });
	turbulenceModels.insert(
			{ "arithmetic3", turbulenceModel::arithmeticA3Model });
	turbulenceModels.insert( { "kEpsA", turbulenceModel::kEpsAModel });
	turbulenceModels.insert( { "BHR", turbulenceModel::BHRModel });
	turbulenceModels.insert( { "BHRKL", turbulenceModel::BHRKLModel });

	turbulenceModel turbulenceModelFlag;
	try
	{
		turbulenceModelFlag = turbulenceModels.at(sourceTypeString);
	} catch (std::out_of_range&)
	{
		throw exception("Unknown source generation flag.",
				errors::initialisationError);
	}

	std::unique_ptr<abstractTurbulenceModel> trbl;
	switch (turbulenceModelFlag)
	{
	case turbulenceModel::decayModel:
		trbl = std::make_unique<decayModel>(meshIn, turbulenceFlag);
		break;
	case turbulenceModel::shearModel:
		trbl = std::make_unique<shearModel>(meshIn, turbulenceFlag);
		break;
	case turbulenceModel::arithmeticA1Model:
	case turbulenceModel::arithmeticA2Model:
	case turbulenceModel::arithmeticA3Model:
		trbl = std::make_unique<arithmeticAModel>(meshIn, turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::kEpsAModel:
		trbl = std::make_unique<kEpsAModel>(meshIn, turbulenceFlag);
		break;
	case turbulenceModel::BHRModel:
		trbl = std::make_unique<BHRModel>(meshIn, parIn, uCellIn, uSurfIn,
				readDataPoint, turbulenceFlag);
		break;
	case turbulenceModel::BHRKLModel:
		trbl = std::make_unique<BHRKLModel>(meshIn, turbulenceFlag);
		break;
	default:
		trbl = std::make_unique<zeroModel>(meshIn, turbulenceFlag);
		break;
	}

	return trbl;
}

void schemi::abstractTurbulenceModel::particlesTimeIntegration(
		const volumeField<vector>&, const surfaceField<vector>&,
		const volumeField<vector>&, const surfaceField<vector>&, const vector&,
		const concentrationsPack<cubicCell>&,
		const std::vector<volumeField<scalar>>&, const boundaryConditionValue&,
		const std::valarray<scalar>&, const scalar, const volumeField<vector>&,
		const volumeField<scalar>&, const volumeField<tensor>&)
{
}

void schemi::abstractTurbulenceModel::particlesWriteOutput(const std::string&,
		const scalar) const
{
}

void schemi::abstractTurbulenceModel::checkTransitionToTurbulenceModel(
		const volumeField<scalar>&, const surfaceField<scalar>&,
		volumeField<scalar>&, volumeField<scalar>&, volumeField<vector>&,
		volumeField<scalar>&, const concentrationsPack<cubicCell>&,
		const boundaryConditionValue&, const scalar) noexcept
{

}

schemi::abstractTurbulenceModel::~abstractTurbulenceModel() noexcept
{
}
