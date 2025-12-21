/*
 * bunchOfFields.hpp
 *
 *  Created on: 2019/12/17
 *      Author: Maxim Boldyrev
 *
 *      Structure storing fields.
 */

#ifndef BUNCHOFFIELDS_HPP_
#define BUNCHOFFIELDS_HPP_

#include <vector>

#include "abstractMixtureThermodynamics.hpp"
#include "mesh.hpp"
#include "scalar.hpp"
#include "vector.hpp"
#include "concentrationsPack.hpp"
#include "field.hpp"
#include "fieldOperations.hpp"

namespace schemi
{
template<typename typeOfEnity1>
struct bunchOfFields
{
	concentrationsPack<typeOfEnity1> concentration;
	std::vector<field<scalar, typeOfEnity1>> density;
	field<vector, typeOfEnity1> velocity;
	field<scalar, typeOfEnity1> pressure;
	field<scalar, typeOfEnity1> kTurb;
	field<scalar, typeOfEnity1> epsTurb;
	field<vector, typeOfEnity1> aTurb;
	field<scalar, typeOfEnity1> bTurb;
	field<vector, typeOfEnity1> momentum;
	field<scalar, typeOfEnity1> internalEnergy;
	field<scalar, typeOfEnity1> totalEnergy;
	field<scalar, typeOfEnity1> temperature;
	field<scalar, typeOfEnity1> rhokTurb;
	field<scalar, typeOfEnity1> rhoepsTurb;
	field<vector, typeOfEnity1> rhoaTurb;
	field<scalar, typeOfEnity1> rhobTurb;
	field<scalar, typeOfEnity1> HelmholtzEnergy;
	field<scalar, typeOfEnity1> entropy;

	bunchOfFields(const mesh & meshRef,
			const std::size_t numberOfcomponents) noexcept :
			concentration { meshRef, numberOfcomponents },

			density { std::size_t(numberOfcomponents + 1), field<scalar,
					typeOfEnity1>(meshRef, 0) },

			velocity { meshRef, vector(0) },

			pressure { meshRef, 0 },

			kTurb { meshRef, scalar { 0 } },

			epsTurb { meshRef, 0 },

			aTurb { meshRef, vector(0) },

			bTurb { meshRef, 0 },

			momentum { meshRef, vector(0) },

			internalEnergy { meshRef, 0 },

			totalEnergy { meshRef, 0 },

			temperature { meshRef, 0 },

			rhokTurb { meshRef, scalar { 0 } },

			rhoepsTurb { meshRef, 0 },

			rhoaTurb { meshRef, vector(0) },

			rhobTurb { meshRef, 0 },

			HelmholtzEnergy { meshRef, 0 },

			entropy { meshRef, 0 }
	{
	}

	template<typename typeOfEnity2>
	explicit bunchOfFields(
			const bunchOfFields<typeOfEnity2> & otherFieldforCopy) noexcept :
			concentration { otherFieldforCopy.pressure.meshRef(),
					otherFieldforCopy.concentration.v.size() - std::size_t(1) },

			density { otherFieldforCopy.density.size(), field<scalar,
					typeOfEnity1> { otherFieldforCopy.pressure.meshRef(), 0 } },

			velocity { otherFieldforCopy.pressure.meshRef(), vector(0),
					otherFieldforCopy.velocity.boundCond() },

			pressure { otherFieldforCopy.pressure.meshRef(), 0,
					otherFieldforCopy.pressure.boundCond() },

			kTurb { otherFieldforCopy.pressure.meshRef(), 0,
					otherFieldforCopy.kTurb.boundCond() },

			epsTurb { otherFieldforCopy.pressure.meshRef(), 0,
					otherFieldforCopy.epsTurb.boundCond() },

			aTurb { otherFieldforCopy.pressure.meshRef(), vector(0),
					otherFieldforCopy.aTurb.boundCond() },

			bTurb { otherFieldforCopy.pressure.meshRef(), 0,
					otherFieldforCopy.bTurb.boundCond() },

			momentum { otherFieldforCopy.pressure.meshRef(), vector(0) },

			internalEnergy { otherFieldforCopy.pressure.meshRef(), 0 },

			totalEnergy { otherFieldforCopy.pressure.meshRef(), 0 },

			temperature { otherFieldforCopy.pressure.meshRef(), 0 },

			rhokTurb { otherFieldforCopy.pressure.meshRef(), 0 },

			rhoepsTurb { otherFieldforCopy.pressure.meshRef(), 0 },

			rhoaTurb { otherFieldforCopy.pressure.meshRef(), vector(0) },

			rhobTurb { otherFieldforCopy.pressure.meshRef(), 0 },

			HelmholtzEnergy { otherFieldforCopy.pressure.meshRef(), 0 },

			entropy { otherFieldforCopy.pressure.meshRef(), 0 }
	{
		for (std::size_t k = 0; k < otherFieldforCopy.concentration.v.size();
				++k)
			if (k != 0)
				concentration.v[k].boundCond_wr() =
						otherFieldforCopy.concentration.v[k].boundCond();
	}

	void clear() noexcept
	{
		for (auto & c : concentration)
			c.v.r() = 0;
		for (auto & d : density)
			d.r() = 0;
		velocity.r() = vector(0);
		pressure.r() = 0;
		kTurb.r() = scalar { 0 };
		epsTurb.r() = 0;
		aTurb.r() = vector(0);
		bTurb.r() = 0;
		momentum.r() = vector(0);
		internalEnergy.r() = 0;
		totalEnergy.r() = 0;
		temperature.r() = 0;
		rhokTurb.r() = scalar { 0 };
		rhoepsTurb.r() = 0;
		rhoaTurb.r() = vector(0);
		rhobTurb.r() = 0;
		HelmholtzEnergy.r() = 0;
		entropy.r() = 0;
	}

	void average(const bunchOfFields<typeOfEnity1> & in,
			const abstractMixtureThermodynamics & mixture,
			const scalar ownWeight) noexcept
	{
		const auto inWeight = 1 - ownWeight;

		for (std::size_t k = 0; k < in.concentration.v.size(); ++k)
		{
			concentration.v[k].val() = (ownWeight * concentration.v[k]
					+ inWeight * in.concentration.v[k]).cval();

			density[k].val() = (ownWeight * density[k]
					+ inWeight * in.density[k]).cval();
		}
		momentum.val() = (momentum * ownWeight + in.momentum * inWeight).cval();
		totalEnergy.val() =
				(ownWeight * totalEnergy + inWeight * in.totalEnergy).cval();
		rhokTurb.val() = (ownWeight * rhokTurb + inWeight * in.rhokTurb).cval();
		rhoepsTurb.val() =
				(ownWeight * rhoepsTurb + inWeight * in.rhoepsTurb).cval();
		rhoaTurb.val() = (rhoaTurb * ownWeight + in.rhoaTurb * inWeight).cval();
		rhobTurb.val() = (ownWeight * rhobTurb + inWeight * in.rhobTurb).cval();

		velocity.val() = (momentum / density[0]).cval();
		kTurb.val() = (rhokTurb / density[0]).cval();
		epsTurb.val() = (rhoepsTurb / density[0]).cval();
		aTurb.val() = (rhoaTurb / density[0]).cval();
		{
			internalEnergy.val() = (totalEnergy - rhokTurb
					- 0.5 * density[0] * (velocity & velocity)).cval();
		}
		pressure.val() = mixture.pFromUv(concentration.p,
				internalEnergy.cval());
		temperature.val() = mixture.TFromUv(concentration.p,
				internalEnergy.cval());
		HelmholtzEnergy.val() = mixture.Fv(concentration.p, temperature.cval());
		entropy.val() = mixture.Sv(concentration.p, temperature.cval());
	}

	void copyFrom(const bunchOfFields<typeOfEnity1> & in,
			const abstractMixtureThermodynamics & mixture) noexcept
	{
		for (std::size_t k = 0; k < in.concentration.v.size(); ++k)
		{
			concentration.v[k].val() = in.concentration.v[k].cval();

			density[k].val() = in.density[k].cval();
		}
		momentum.val() = in.momentum.cval();
		totalEnergy.val() = in.totalEnergy.cval();
		rhokTurb.val() = in.rhokTurb.cval();
		rhoepsTurb.val() = in.rhoepsTurb.cval();
		rhoaTurb.val() = in.rhoaTurb.cval();
		rhobTurb.val() = in.rhobTurb.cval();

		velocity.val() = (momentum / density[0]).cval();
		kTurb.val() = (rhokTurb / density[0]).cval();
		epsTurb.val() = (rhoepsTurb / density[0]).cval();
		aTurb.val() = (rhoaTurb / density[0]).cval();
		{
			internalEnergy.val() = (totalEnergy - rhokTurb
					- 0.5 * density[0] * (velocity & velocity)).cval();
		}
		pressure.val() = mixture.pFromUv(concentration.p,
				internalEnergy.cval());
		temperature.val() = mixture.TFromUv(concentration.p,
				internalEnergy.cval());
		HelmholtzEnergy.val() = mixture.Fv(concentration.p, temperature.cval());
		entropy.val() = mixture.Sv(concentration.p, temperature.cval());
	}
};
}  // namespace schemi

#endif /* BUNCHOFFIELDS_HPP_ */
