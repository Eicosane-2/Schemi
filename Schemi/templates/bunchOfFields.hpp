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
#include "fieldProducts.hpp"

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
			concentration.v[k].r() = ownWeight * concentration.v[k]()
					+ inWeight * in.concentration.v[k]();

			density[k].r() = ownWeight * density[k]()
					+ inWeight * in.density[k]();
		}
		momentum.r() = astProduct(momentum, ownWeight)()
				+ astProduct(in.momentum, inWeight)();
		totalEnergy.r() = ownWeight * totalEnergy()
				+ inWeight * in.totalEnergy();
		rhokTurb.r() = ownWeight * rhokTurb() + inWeight * in.rhokTurb();
		rhoepsTurb.r() = ownWeight * rhoepsTurb() + inWeight * in.rhoepsTurb();
		rhoaTurb.r() = astProduct(rhoaTurb, ownWeight)()
				+ astProduct(in.rhoaTurb, inWeight)();
		rhobTurb.r() = ownWeight * rhobTurb() + inWeight * in.rhobTurb();

		velocity.r() = division(momentum, density[0])();
		kTurb.r() = rhokTurb() / density[0]();
		epsTurb.r() = rhoepsTurb() / density[0]();
		aTurb.r() = division(rhoaTurb, density[0])();
		{
			internalEnergy.r() = totalEnergy() - rhokTurb()
					- 0.5 * density[0]() * ampProduct(velocity, velocity)();
		}
		pressure.r() = mixture.pFromUv(concentration.p, internalEnergy());
		temperature.r() = mixture.TFromUv(concentration.p, internalEnergy());
		HelmholtzEnergy.r() = mixture.Fv(concentration.p, temperature());
		entropy.r() = mixture.Sv(concentration.p, temperature());
	}

	void copyFrom(const bunchOfFields<typeOfEnity1> & in,
			const abstractMixtureThermodynamics & mixture) noexcept
	{
		for (std::size_t k = 0; k < in.concentration.v.size(); ++k)
		{
			concentration.v[k].r() = in.concentration.v[k]();

			density[k].r() = in.density[k]();
		}
		momentum.r() = in.momentum();
		totalEnergy.r() = in.totalEnergy();
		rhokTurb.r() = in.rhokTurb();
		rhoepsTurb.r() = in.rhoepsTurb();
		rhoaTurb.r() = in.rhoaTurb();
		rhobTurb.r() = in.rhobTurb();

		velocity.r() = division(momentum, density[0])();
		kTurb.r() = rhokTurb() / density[0]();
		epsTurb.r() = rhoepsTurb() / density[0]();
		aTurb.r() = division(rhoaTurb, density[0])();
		{
			internalEnergy.r() = totalEnergy() - rhokTurb()
					- 0.5 * density[0]() * ampProduct(velocity, velocity)();
		}
		pressure.r() = mixture.pFromUv(concentration.p, internalEnergy());
		temperature.r() = mixture.TFromUv(concentration.p, internalEnergy());
		HelmholtzEnergy.r() = mixture.Fv(concentration.p, temperature());
		entropy.r() = mixture.Sv(concentration.p, temperature());
	}
};
}  // namespace schemi

#endif /* BUNCHOFFIELDS_HPP_ */
