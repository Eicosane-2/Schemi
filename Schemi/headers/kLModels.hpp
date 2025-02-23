/*
 * kLModels.hpp
 *
 *  Created on: 2024/12/03
 *      Author: Maxim Boldyrev
 */

#ifndef KLMODELS_HPP_
#define KLMODELS_HPP_

#include "abstractTurbulenceModel.hpp"

namespace schemi
{
class kLModels: public abstractTurbulenceModel
{
	scalar CMS_R_;
	scalar CMS_D_;

	scalar thetaS(const scalar divV, const scalar k, const scalar L,
			const scalar CMS) const noexcept;

protected:
	void CMS_RSet(const scalar) noexcept;
	void CMS_DSet(const scalar) noexcept;

public:
	virtual ~kLModels() noexcept =0;

	scalar CMS_R() const noexcept;
	scalar CMS_D() const noexcept;

	explicit kLModels(const mesh & meshIn, const bool turb = true,
			const bool a = false, const bool b = false,
			const turbulenceModel model = turbulenceModel::unknownModel,
			const scalar Cmu_in = 0.28, const scalar sigmaSc_in = 0.5,
			const scalar sigmaT_in = 1, const scalar sigmaE_in = 1,
			const scalar sigmak_in = 0.5, const scalar sigmaL_in = 0.1,
			const scalar sigmaa_in = 0.5, const scalar sigmab_in = 0.5,
			const scalar CMS_R_in = 0.1, const scalar CMS_D_in = 1.0) noexcept;

	std::valarray<scalar> calculateNut(const std::valarray<scalar> & k,
			const std::valarray<scalar> & L) const noexcept override;

	scalar calculateNut(const scalar k, const scalar L) const noexcept override;

	virtual std::valarray<scalar> rhoepsilon(
			const bunchOfFields<cubicCell> & cf,
			const abstractMixtureThermodynamics & th) const noexcept override;

	scalar thetaS_R(const scalar divV, const scalar k,
			const scalar L) const noexcept override;

	volumeField<scalar> thetaS_R(const volumeField<scalar> & divV,
			const volumeField<scalar> & k,
			const volumeField<scalar> & L) const noexcept override;

	surfaceField<scalar> thetaS_R(const surfaceField<scalar> & divV,
			const surfaceField<scalar> & k,
			const surfaceField<scalar> & L) const noexcept override;

	scalar thetaS_D(const scalar divV, const scalar k,
			const scalar L) const noexcept override;

	volumeField<scalar> thetaS_D(const volumeField<scalar> & divV,
			const volumeField<scalar> & k,
			const volumeField<scalar> & L) const noexcept override;

	surfaceField<scalar> thetaS_D(const surfaceField<scalar> & divV,
			const surfaceField<scalar> & k,
			const surfaceField<scalar> & L) const noexcept override;
};
}  // namespace schemi

#endif /* KLMODELS_HPP_ */
