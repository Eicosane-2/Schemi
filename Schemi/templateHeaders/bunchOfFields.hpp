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

#include "concentrationsPack.hpp"
#include "field.hpp"
#include "fieldProducts.hpp"
#include "abstractMixtureThermodynamics.hpp"
#include "mesh.hpp"
#include "scalar.hpp"
#include "vector.hpp"

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
	bunchOfFields(
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
				concentration.v[k].boundCond_r() =
						otherFieldforCopy.concentration.v[k].boundCond();
	}

	void clear() noexcept
	{
		for (auto & c : concentration)
			c.v.ref_r() = 0;
		for (auto & d : density)
			d.ref_r() = 0;
		velocity.ref_r() = vector(0);
		pressure.ref_r() = 0;
		kTurb.ref_r() = scalar { 0 };
		epsTurb.ref_r() = 0;
		aTurb.ref_r() = vector(0);
		bTurb.ref_r() = 0;
		momentum.ref_r() = vector(0);
		internalEnergy.ref_r() = 0;
		totalEnergy.ref_r() = 0;
		temperature.ref_r() = 0;
		rhokTurb.ref_r() = scalar { 0 };
		rhoepsTurb.ref_r() = 0;
		rhoaTurb.ref_r() = vector(0);
		rhobTurb.ref_r() = 0;
		HelmholtzEnergy.ref_r() = 0;
		entropy.ref_r() = 0;
	}

	void average(const bunchOfFields<typeOfEnity1> & in,
			const abstractMixtureThermodynamics & mixture,
			const scalar ownWeight) noexcept
	{
		const auto inWeight = 1 - ownWeight;

		for (std::size_t k = 0; k < in.concentration.v.size(); ++k)
		{
			concentration.v[k].ref_r() = ownWeight * concentration.v[k].ref()
					+ inWeight * in.concentration.v[k].ref();

			density[k].ref_r() = ownWeight * density[k].ref()
					+ inWeight * in.density[k].ref();
		}
		momentum.ref_r() = astProduct(momentum, ownWeight).ref()
				+ astProduct(in.momentum, inWeight).ref();
		totalEnergy.ref_r() = ownWeight * totalEnergy.ref()
				+ inWeight * in.totalEnergy.ref();
		rhokTurb.ref_r() = ownWeight * rhokTurb.ref()
				+ inWeight * in.rhokTurb.ref();
		rhoepsTurb.ref_r() = ownWeight * rhoepsTurb.ref()
				+ inWeight * in.rhoepsTurb.ref();
		rhoaTurb.ref_r() = astProduct(rhoaTurb, ownWeight).ref()
				+ astProduct(in.rhoaTurb, inWeight).ref();
		rhobTurb.ref_r() = ownWeight * rhobTurb.ref()
				+ inWeight * in.rhobTurb.ref();

		velocity.ref_r() = division(momentum, density[0]).ref();
		kTurb.ref_r() = rhokTurb.ref() / density[0].ref();
		epsTurb.ref_r() = rhoepsTurb.ref() / density[0].ref();
		aTurb.ref_r() = division(rhoaTurb, density[0]).ref();
		{
			internalEnergy.ref_r() = totalEnergy.ref() - rhokTurb.ref()
					- 0.5 * density[0].ref()
							* ampProduct(velocity, velocity).ref();
		}
		pressure.ref_r() = mixture.pFromUv(concentration.p,
				internalEnergy.ref());
		temperature.ref_r() = mixture.TFromUv(concentration.p,
				internalEnergy.ref());
		HelmholtzEnergy.ref_r() = mixture.Fv(concentration.p,
				temperature.ref());
		entropy.ref_r() = mixture.Sv(concentration.p, temperature.ref());
	}

	void copyFrom(const bunchOfFields<typeOfEnity1> & in,
			const abstractMixtureThermodynamics & mixture) noexcept
	{
		for (std::size_t k = 0; k < in.concentration.v.size(); ++k)
		{
			concentration.v[k].ref_r() = in.concentration.v[k].ref();

			density[k].ref_r() = in.density[k].ref();
		}
		momentum.ref_r() = in.momentum.ref();
		totalEnergy.ref_r() = in.totalEnergy.ref();
		rhokTurb.ref_r() = in.rhokTurb.ref();
		rhoepsTurb.ref_r() = in.rhoepsTurb.ref();
		rhoaTurb.ref_r() = in.rhoaTurb.ref();
		rhobTurb.ref_r() = in.rhobTurb.ref();

		velocity.ref_r() = division(momentum, density[0]).ref();
		kTurb.ref_r() = rhokTurb.ref() / density[0].ref();
		epsTurb.ref_r() = rhoepsTurb.ref() / density[0].ref();
		aTurb.ref_r() = division(rhoaTurb, density[0]).ref();
		{
			internalEnergy.ref_r() = totalEnergy.ref() - rhokTurb.ref()
					- 0.5 * density[0].ref()
							* ampProduct(velocity, velocity).ref();
		}
		pressure.ref_r() = mixture.pFromUv(concentration.p,
				internalEnergy.ref());
		temperature.ref_r() = mixture.TFromUv(concentration.p,
				internalEnergy.ref());
		HelmholtzEnergy.ref_r() = mixture.Fv(concentration.p,
				temperature.ref());
		entropy.ref_r() = mixture.Sv(concentration.p, temperature.ref());
	}
};
}  // namespace schemi

#endif /* BUNCHOFFIELDS_HPP_ */
