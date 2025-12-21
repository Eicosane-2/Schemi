/*
 * abstractTurbulenceModel.hpp
 *
 *  Created on: 2024/12/03
 *      Author: Maxim Boldyrev
 */

#ifndef ABSTRACTTURBULENCEMODEL_HPP_
#define ABSTRACTTURBULENCEMODEL_HPP_

#include <memory>

#include "bunchOfFields.hpp"
#include "diffusiveFields.hpp"
#include "volumeField.hpp"
#include "surfaceField.hpp"
#include "turbulenceModelEnum.hpp"
#include "MPIHandler.hpp"

namespace schemi
{
class boundaryConditionValue;

class abstractTurbulenceModel
{
	bool turbulence_;
	const bool aField_, bField_;
	const turbulenceModel modelType_;

	scalar Cmu_;
	scalar sigmaSc_;
	scalar sigmaT_;
	scalar sigmaE_;
	scalar sigmaK_;
	scalar sigmaEps_;
	scalar sigmaa_;
	scalar sigmab_;

	std::valarray<scalar> modelParameters_;

	volumeField<scalar> y;

	std::valarray<scalar> solveMatrix(const volumeField<scalar> & phi,
			const boundaryConditionValue & boundCond) const;

	scalar mink_value_ { 1E-15 };
	scalar minepsilon_value_ { 1E-15 };
	scalar minb_value_ { 1E-15 };

protected:
	void turbulenceSet(const bool) noexcept;
	void CmuSet(const scalar) noexcept;
	void sigmaScSet(const scalar) noexcept;
	void sigmaTSet(const scalar) noexcept;
	void sigmaESet(const scalar) noexcept;
	void sigmaKSet(const scalar) noexcept;
	void sigmaEpsSet(const scalar) noexcept;
	void sigmaaSet(const scalar) noexcept;
	void sigmabSet(const scalar) noexcept;

	std::valarray<scalar>& modelParametersSet() noexcept;

	const std::valarray<scalar>& modelParameters() const noexcept;

	void minkSet(const scalar) noexcept;
	void minepsilonSet(const scalar) noexcept;
	void minbSet(const scalar) noexcept;

	constexpr static scalar aFlag { 1E15 };

public:
	virtual ~abstractTurbulenceModel() noexcept =0;

	bool turbulence() const noexcept;
	bool aField() const noexcept;
	bool bField() const noexcept;
	turbulenceModel model() const noexcept;
	scalar Cmu() const noexcept;
	scalar sigmaSc() const noexcept;
	scalar sigmaT() const noexcept;
	scalar sigmaE() const noexcept;
	scalar sigmaK() const noexcept;
	scalar sigmaEps() const noexcept;
	scalar sigmaa() const noexcept;
	scalar sigmab() const noexcept;
	scalar mink() const noexcept;
	scalar minepsilon() const noexcept;
	scalar minb() const noexcept;

	abstractTurbulenceModel(const mesh & meshIn, const bool turb = false,
			const bool a = false, const bool b = false,
			const turbulenceModel model = turbulenceModel::unknownModel,
			const scalar Cmu_in = 0.1, const scalar sigmaSc_in = 1,
			const scalar sigmaT_in = 1, const scalar sigmaE_in = 1,
			const scalar sigmaK_in = 1, const scalar sigmaEps_in = 1,
			const scalar sigmaa_in = 1, const scalar sigmab_in = 1) noexcept;

	virtual std::valarray<scalar> calculateNut(
			const std::valarray<scalar>& /*k*/,
			const std::valarray<scalar>& /*eps*/) const noexcept =0;

	virtual scalar calculateNut(const scalar /*k*/,
			const scalar /*eps*/) const noexcept =0;

	virtual std::valarray<scalar> rhoepsilon(
			const bunchOfFields<cubicCell>& /*cf*/,
			const abstractMixtureThermodynamics& /*th*/,
			const volumeField<scalar>& /*k*/,
			const volumeField<scalar>& /*eps*/) const noexcept =0;

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

	virtual std::tuple<std::pair<volumeField<scalar>, volumeField<scalar>>,
			std::pair<volumeField<scalar>, volumeField<scalar>>,
			std::pair<volumeField<vector>, volumeField<vector>>,
			std::pair<volumeField<scalar>, volumeField<scalar>>,
			volumeField<scalar>> calculate(scalar & sourceTimestep,
			const scalar sourceTimestepCoeff,
			const bunchOfFields<cubicCell> & cellFields,
			const diffusiveFields & diffFieldsOld,
			const volumeField<tensor> & gradV,
			const volumeField<vector> & divDevPhysVisc,
			const volumeField<vector> & gradP,
			const volumeField<vector> & gradRho,
			const volumeField<tensor> & grada, const volumeField<scalar> & diva,
			const volumeField<vector> & gradb,
			const volumeField<tensor> & spherR,
			const volumeField<tensor> & devR,
			const volumeField<vector> & gradMav_n,
			const abstractMixtureThermodynamics & mixture,
			const volumeField<scalar> & nu_t) const noexcept =0;

	void calculateNearWallDistance(const volumeField<scalar> & eps,
			const boundaryConditionValue & boundCond);

	static std::unique_ptr<abstractTurbulenceModel> createTurbulenceModel(
			const mesh & meshIn, const std::string & turbulenceONString,
			const std::string & sourceTypeString, const MPIHandler & parIn,
			const volumeField<vector> & uCellIn,
			const surfaceField<vector> & uSurfIn,
			const std::pair<std::size_t, std::string> & readDataPoint);

	virtual void particlesTimeIntegration(
			const volumeField<vector> & gradRhoCell,
			const surfaceField<vector> & gradRhoSurf,
			const volumeField<vector> & uCell,
			const surfaceField<vector> & uSurf,
			const concentrationsPack<cubicCell> & concentrations,
			const std::vector<volumeField<scalar>> & densities,
			const boundaryConditionValue & boundVal,
			const std::valarray<scalar> & M, const scalar timestep,
			const volumeField<vector> & gradP, const volumeField<scalar> & divU,
			const volumeField<tensor> & gradU);

	virtual void particlesWriteOutput(
			const std::string & fieldDataDirectoryName,
			const scalar Time) const;

	virtual void checkTransitionToTurbulenceModel(
			const volumeField<scalar> & nuCell,
			const surfaceField<scalar> & nuSurface, volumeField<scalar> & k,
			volumeField<scalar> & epsilon, volumeField<vector> & a,
			volumeField<scalar> & b,
			const concentrationsPack<cubicCell> & concentrations,
			const boundaryConditionValue & boundVal,
			const scalar timestep) noexcept;

	virtual bool isInitialisationModelUsed() const noexcept
	{
		return false;
	}
};
}  // namespace schemi

#endif /* ABSTRACTTURBULENCEMODEL_HPP_ */
