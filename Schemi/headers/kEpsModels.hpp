/*
 * kEpsModels.hpp
 *
 *  Created on: 2024/12/03
 *      Author: Maxim Boldyrev
 */

#ifndef KEPSMODELS_HPP_
#define KEPSMODELS_HPP_

#include "abstractTurbulenceModel.hpp"

namespace schemi
{
class kEpsModels: public abstractTurbulenceModel
{
	scalar CMS_R_;
	scalar CMS_D_;

	scalar thetaS(const scalar divV, const scalar k, const scalar eps,
			const scalar CMS) const noexcept;

protected:
	void CMS_RSet(const scalar) noexcept;
	void CMS_DSet(const scalar) noexcept;

public:
	virtual ~kEpsModels() noexcept =0;

	scalar CMS_R() const noexcept;
	scalar CMS_D() const noexcept;

	explicit kEpsModels(const mesh & meshIn, const bool turb = true,
			const bool a = false, const bool b = false,
			const turbulenceModel model = turbulenceModel::unknownModel,
			const scalar Cmu_in = 0.09, const scalar sigmaSc_in = 0.5,
			const scalar sigmaT_in = 1, const scalar sigmaE_in = 1,
			const scalar sigmak_in = 0.5, const scalar sigmaeps_in = 0.4,
			const scalar sigmaa_in = 0.5, const scalar sigmab_in = 0.5,
			const scalar CMS_R_in = 0.1, const scalar CMS_D_in = 1.0) noexcept;

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

#endif /* KEPSMODELS_HPP_ */
