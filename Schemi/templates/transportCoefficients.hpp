/*
 * transportCoefficients.hpp
 *
 *  Created on: 2020/01/25
 *      Author: Maxim Boldyrev
 *
 *      Structure storing fields for diffusion coefficients.
 */

#ifndef TRANSPORTCOEFFICIENTS_HPP_
#define TRANSPORTCOEFFICIENTS_HPP_

#include "turbulenceModelEnum.hpp"
#include "abstractTurbulentParameters.hpp"

namespace schemi
{
template<typename typeOfEntity>
struct transportCoefficients
{
	field<scalar, typeOfEntity> tNu;
	field<scalar, typeOfEntity> tD;
	field<scalar, typeOfEntity> tKappa;
	field<scalar, typeOfEntity> tLambda;
	field<scalar, typeOfEntity> k_D;
	field<scalar, typeOfEntity> eps_D;
	field<scalar, typeOfEntity> a_D;
	field<scalar, typeOfEntity> b_D;

	field<scalar, typeOfEntity> pNu;
	field<scalar, typeOfEntity> pD;
	field<scalar, typeOfEntity> pKappa;

	explicit transportCoefficients(const mesh & meshRef) noexcept :
			tNu { meshRef, 0 },

			tD { meshRef, 0 },

			tKappa { meshRef, 0 },

			tLambda { meshRef, 0 },

			k_D { meshRef, 0 },

			eps_D { meshRef, 0 },

			a_D { meshRef, 0 },

			b_D { meshRef, 0 },

			pNu { meshRef, 0 },

			pD { meshRef, 0 },

			pKappa { meshRef, 0 }
	{
	}

	void recalculateCoefficients(const field<scalar, typeOfEntity> & kField,
			const field<scalar, typeOfEntity> & epsField,
			const abstractTurbulentParameters & tp) noexcept
	{
		tNu.ref_r() = tp.calculateNut(kField.ref(), epsField.ref());

		recalculateCoefficients(tp);
	}

	void recalculateCoefficients(
			const abstractTurbulentParameters & tp) noexcept
	{
		tD.ref_r() = tNu.ref() / tp.sigmaSc();
		tKappa.ref_r() = tNu.ref() / tp.sigmaT();
		tLambda.ref_r() = tNu.ref() / tp.sigmaE();
		k_D.ref_r() = tNu.ref() / tp.sigmaK();
		eps_D.ref_r() = tNu.ref() / tp.sigmaEps();
		a_D.ref_r() = tNu.ref() / tp.sigmaa();
		b_D.ref_r() = tNu.ref() / tp.sigmab();
	}

	void tAssign(const scalar val) noexcept
	{
		tNu.ref_r() = val;
		tD.ref_r() = val;
		tKappa.ref_r() = val;
		tLambda.ref_r() = val;
		k_D.ref_r() = val;
		eps_D.ref_r() = val;
		a_D.ref_r() = val;
		b_D.ref_r() = val;
	}

	auto& operator=(const scalar val) noexcept
	{
		tNu.ref_r() = val;
		tD.ref_r() = val;
		tKappa.ref_r() = val;
		tLambda.ref_r() = val;
		k_D.ref_r() = val;
		eps_D.ref_r() = val;
		a_D.ref_r() = val;
		b_D.ref_r() = val;

		pNu.ref_r() = val;
		pD.ref_r() = val;
		pKappa.ref_r() = val;

		return *this;
	}
};

template<typename typeOfEnity>
struct effectiveTransportCoefficients: transportCoefficients<typeOfEnity>
{
	field<scalar, typeOfEnity> mu;
	std::vector<field<scalar, typeOfEnity>> mD;
	field<scalar, typeOfEnity> rhoD;
	field<scalar, typeOfEnity> kappa;
	field<scalar, typeOfEnity> rhoDk;
	field<scalar, typeOfEnity> rhoDeps;
	field<scalar, typeOfEnity> rhoDa;
	field<scalar, typeOfEnity> rhoDb;

	effectiveTransportCoefficients(const mesh & meshRef,
			const std::size_t compt) noexcept :
			transportCoefficients<typeOfEnity>(meshRef), mu { meshRef, 0 },

			mD { compt, field<scalar, typeOfEnity>(meshRef, 0) },

			rhoD { meshRef, 0 },

			kappa { meshRef, 0 },

			rhoDk { meshRef, 0 },

			rhoDeps { meshRef, 0 },

			rhoDa { meshRef, 0 },

			rhoDb { meshRef, 0 }
	{
	}

	field<scalar, typeOfEnity> maxValue(const bool turbulenceFlag,
			const turbulenceModelEnum sourceFlag) const noexcept
	{
		auto maxFieldVal = this->pNu;

		for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
			maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
					this->pD.ref()[i]);

		if (turbulenceFlag)
		{
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						this->tNu.ref()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						this->tD.ref()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						2 * this->tLambda.ref()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						this->k_D.ref()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						this->eps_D.ref()[i]);

			if ((sourceFlag == turbulenceModelEnum::BHRSource)
					|| (sourceFlag == turbulenceModelEnum::BHRKLSource)
					|| (sourceFlag == turbulenceModelEnum::kEpsASource))
			{
				for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
					maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
							this->a_D.ref()[i]);

				if ((sourceFlag == turbulenceModelEnum::BHRSource)
						|| (sourceFlag == turbulenceModelEnum::BHRKLSource))
					for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
						maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
								this->b_D.ref()[i]);
			}
		}

		return maxFieldVal;
	}
};
}  // namespace schemi

#endif /* TRANSPORTCOEFFICIENTS_HPP_ */
