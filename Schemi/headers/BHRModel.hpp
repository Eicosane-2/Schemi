/*
 * BHRModel.hpp
 *
 *  Created on: 2025/03/25
 *      Author: Maxim Boldyrev
 */

#ifndef BHRMODEL_HPP_
#define BHRMODEL_HPP_

#include "kEpsModels.hpp"
#include "diffusiveFields.hpp"
#include "instabilityParticlesHandler.hpp"

namespace schemi
{
class BHRModel: public kEpsModels
{
	scalar C0() const noexcept;
	scalar C1() const noexcept;
	scalar C2() const noexcept;
	scalar C3() const noexcept;
	scalar C4() const noexcept;

	scalar Ca() const noexcept;
	scalar Cb1() const noexcept;
	scalar Cb2() const noexcept;

	scalar alpha2() const noexcept;
	scalar alpha3() const noexcept;
	scalar alpha4() const noexcept;

	scalar CMSM() const noexcept;

	instabilityParticlesHandler initialisation;

public:
	BHRModel(const mesh & meshIn, const MPIHandler & parIn,
			const volumeField<vector> & uCellIn,
			const surfaceField<vector> & uSurfIn,
			const std::pair<std::size_t, std::string> & readDataPoint,
			const bool turb_in);

	std::tuple<std::pair<volumeField<scalar>, volumeField<scalar>>,
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
			const volumeField<scalar> & nu_t) const noexcept override;

	std::valarray<scalar> rhoepsilon(const bunchOfFields<cubicCell> & cf,
			const abstractMixtureThermodynamics & th,
			const volumeField<scalar> & k,
			const volumeField<scalar> & eps) const noexcept override;

	void particlesTimeIntegration(const volumeField<vector> & gradRhoCell,
			const surfaceField<vector> & gradRhoSurf,
			const volumeField<vector> & uCell,
			const surfaceField<vector> & uSurf, const vector & g,
			const concentrationsPack<cubicCell> & concentrations,
			const std::vector<volumeField<scalar>> & densities,
			const boundaryConditionValue & boundVal,
			const std::valarray<scalar> & M, const scalar timestep,
			const volumeField<vector> & gradP, const volumeField<scalar> & divU,
			const volumeField<tensor> & gradU) override;

	void particlesWriteOutput(const std::string & fieldDataDirectoryName,
			const scalar Time) const override;

	void checkTransitionToTurbulenceModel(const volumeField<scalar> & nuCell,
			const surfaceField<scalar> & nuSurface, volumeField<scalar> & k,
			volumeField<scalar> & epsilon, volumeField<vector> & a,
			volumeField<scalar> & b,
			const concentrationsPack<cubicCell> & concentrations,
			const boundaryConditionValue & boundVal,
			const scalar timestep) noexcept override;

	bool isInitialisationModelUsed() const noexcept override
	{
		return initialisation.isInitialisationModelUsed();
	}
};
}  // namespace schemi

#endif /* BHRMODEL_HPP_ */
