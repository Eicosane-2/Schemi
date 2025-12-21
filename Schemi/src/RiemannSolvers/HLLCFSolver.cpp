/*
 * HLLCFSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HLLCFSolver.hpp"

#include "vector.hpp"

std::tuple<schemi::conservativeFlows, schemi::starFields> schemi::HLLCFSolver::calculateFlows(
		const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
		const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
{
	auto & mesh_ { surfaceOwnerSide.pressure.meshRef() };

	conservativeFlows numFluxes(mesh_,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());
	starFields starValues(mesh_,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());

	errors errValue { errors::noErrors };

	const std::valarray<scalar> sonicSpeedOwner(
			std::sqrt(
					surfaceOwnerSide.phaseThermodynamics->sqSonicSpeed(
							surfaceOwnerSide.concentration.p,
							surfaceOwnerSide.density[0].cval(),
							surfaceOwnerSide.internalEnergy.cval(),
							surfaceOwnerSide.pressure.cval())
							+ twothirds * surfaceOwnerSide.kTurb.cval()
									* (scalar(1.)
											+ surfaceOwnerSide.phaseThermodynamics->dpdUv(
													surfaceOwnerSide.concentration.p,
													surfaceOwnerSide.internalEnergy.cval()))));

	const std::valarray<scalar> sonicSpeedNeighbour(
			std::sqrt(
					surfaceNeighbourSide.phaseThermodynamics->sqSonicSpeed(
							surfaceNeighbourSide.concentration.p,
							surfaceNeighbourSide.density[0].cval(),
							surfaceNeighbourSide.internalEnergy.cval(),
							surfaceNeighbourSide.pressure.cval())
							+ twothirds * surfaceNeighbourSide.kTurb.cval()
									* (scalar(1.)
											+ surfaceNeighbourSide.phaseThermodynamics->dpdUv(
													surfaceNeighbourSide.concentration.p,
													surfaceNeighbourSide.internalEnergy.cval()))));

	for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
	{
		const scalar velocityProjectionOwner {
				surfaceOwnerSide.velocity.cval()[i] & mesh_.surfaces()[i].N() };
		const scalar velocityProjectionNeighbour {
				surfaceNeighbourSide.velocity.cval()[i]
						& mesh_.surfaces()[i].N() };

		const scalar turbulentPressureOwner { twothirds
				* surfaceOwnerSide.rhokTurb.cval()[i] };
		const scalar turbulentPressureNeighbour { twothirds
				* surfaceNeighbourSide.rhokTurb.cval()[i] };

		const scalar SOwner { std::min(
				velocityProjectionOwner - sonicSpeedOwner[i],
				velocityProjectionNeighbour - sonicSpeedNeighbour[i]) };
		const scalar SNeighbour { std::max(
				velocityProjectionOwner + sonicSpeedOwner[i],
				velocityProjectionNeighbour + sonicSpeedNeighbour[i]) };

		const scalar Sc { ((surfaceNeighbourSide.pressure.cval()[i]
				+ turbulentPressureNeighbour)
				- (surfaceOwnerSide.pressure.cval()[i] + turbulentPressureOwner)
				+ surfaceOwnerSide.density[0].cval()[i]
						* velocityProjectionOwner
						* (SOwner - velocityProjectionOwner)
				- surfaceNeighbourSide.density[0].cval()[i]
						* velocityProjectionNeighbour
						* (SNeighbour - velocityProjectionNeighbour))
				/ (surfaceOwnerSide.density[0].cval()[i]
						* (SOwner - velocityProjectionOwner)
						- surfaceNeighbourSide.density[0].cval()[i]
								* (SNeighbour - velocityProjectionNeighbour)) };

		const scalar POwnNei {
				0.5
						* ((surfaceOwnerSide.pressure.cval()[i]
								+ turbulentPressureOwner)
								+ (surfaceNeighbourSide.pressure.cval()[i]
										+ turbulentPressureNeighbour)
								+ surfaceOwnerSide.density[0].cval()[i]
										* (SOwner - velocityProjectionOwner)
										* (Sc - velocityProjectionOwner)
								+ surfaceNeighbourSide.density[0].cval()[i]
										* (SNeighbour
												- velocityProjectionNeighbour)
										* (Sc - velocityProjectionNeighbour)) };

		std::valarray<scalar> densityState(numFluxes.density.size());
		vector momentumState;
		scalar totalEnergyState;
		scalar rhokState;
		scalar rhoEpsState;
		vector rhoaState;
		scalar rhobState;
		vector velocityState;
		scalar pressureState;
		vector aState;
		scalar bState;

		if (SOwner >= 0)
		{
			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceOwnerSide.density[k].cval()[i];

			momentumState = surfaceOwnerSide.momentum.cval()[i];
			totalEnergyState = surfaceOwnerSide.totalEnergy.cval()[i];
			rhokState = surfaceOwnerSide.rhokTurb.cval()[i];
			rhoEpsState = surfaceOwnerSide.rhoepsTurb.cval()[i];
			rhoaState = surfaceOwnerSide.rhoaTurb.cval()[i];
			rhobState = surfaceOwnerSide.rhobTurb.cval()[i];

			velocityState = surfaceOwnerSide.velocity.cval()[i];
			pressureState = surfaceOwnerSide.pressure.cval()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];
		}
		else if (Sc >= 0)
		{
			const scalar factor { SOwner - velocityProjectionOwner };
			const scalar denominator { 1 / (SOwner - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceOwnerSide.density[k].cval()[i] * factor
						* denominator;

			momentumState = (surfaceOwnerSide.momentum.cval()[i] * factor
					- mesh_.surfaces()[i].N()
							* (surfaceOwnerSide.pressure.cval()[i]
									+ turbulentPressureOwner - POwnNei))
					* denominator;
			totalEnergyState = (surfaceOwnerSide.totalEnergy.cval()[i] * factor
					- (surfaceOwnerSide.pressure.cval()[i]
							+ turbulentPressureOwner) * velocityProjectionOwner
					+ POwnNei * Sc) * denominator;
			rhokState = surfaceOwnerSide.rhokTurb.cval()[i] * factor
					* denominator;
			rhoEpsState = surfaceOwnerSide.rhoepsTurb.cval()[i] * factor
					* denominator;
			rhoaState = surfaceOwnerSide.rhoaTurb.cval()[i] * factor
					* denominator;
			rhobState = surfaceOwnerSide.rhobTurb.cval()[i] * factor
					* denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceOwnerSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];
		}
		else if (SNeighbour > 0)
		{
			const scalar factor { SNeighbour - velocityProjectionNeighbour };
			const scalar denominator { 1 / (SNeighbour - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceNeighbourSide.density[k].cval()[i]
						* factor * denominator;

			momentumState = (surfaceNeighbourSide.momentum.cval()[i] * factor
					- mesh_.surfaces()[i].N()
							* (surfaceNeighbourSide.pressure.cval()[i]
									+ turbulentPressureNeighbour - POwnNei))
					* denominator;
			totalEnergyState = (surfaceNeighbourSide.totalEnergy.cval()[i]
					* factor
					- (surfaceNeighbourSide.pressure.cval()[i]
							+ turbulentPressureNeighbour)
							* velocityProjectionNeighbour + POwnNei * Sc)
					* denominator;
			rhokState = surfaceNeighbourSide.rhokTurb.cval()[i] * factor
					* denominator;
			rhoEpsState = surfaceNeighbourSide.rhoepsTurb.cval()[i] * factor
					* denominator;
			rhoaState = surfaceNeighbourSide.rhoaTurb.cval()[i] * factor
					* denominator;
			rhobState = surfaceNeighbourSide.rhobTurb.cval()[i] * factor
					* denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceOwnerSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];
		}
		else if (SNeighbour <= 0)
		{
			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceNeighbourSide.density[k].cval()[i];

			momentumState = surfaceNeighbourSide.momentum.cval()[i];
			totalEnergyState = surfaceNeighbourSide.totalEnergy.cval()[i];
			rhokState = surfaceNeighbourSide.rhokTurb.cval()[i];
			rhoEpsState = surfaceNeighbourSide.rhoepsTurb.cval()[i];
			rhoaState = surfaceNeighbourSide.rhoaTurb.cval()[i];
			rhobState = surfaceNeighbourSide.rhobTurb.cval()[i];

			velocityState = surfaceNeighbourSide.velocity.cval()[i];
			pressureState = surfaceNeighbourSide.pressure.cval()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];
		}
		else
			[[unlikely]]
			throw exception("Fatal error in Riemann solver.", errValue =
					errors::RiemannSolverError);

		for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
			numFluxes.density[k].val()[i] = velocityState * densityState[k];

		numFluxes.momentum.val()[i] = momentumState * velocityState
				+ tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
						* (pressureState + twothirds * rhokState);

		numFluxes.totalEnergy.val()[i] = velocityState * totalEnergyState
				+ velocityState * (pressureState + twothirds * rhokState);

		numFluxes.rhokTurb.val()[i] = rhokState * velocityState;
		numFluxes.rhoepsTurb.val()[i] = velocityState * rhoEpsState;
		numFluxes.rhoaTurb.val()[i] = rhoaState * velocityState;
		numFluxes.rhobTurb.val()[i] = velocityState * rhobState;

		starValues.c.v[0].val()[i] = 0;
		for (std::size_t k = 1; k < densityState.size(); ++k)
		{
			starValues.c.v[k].val()[i] = densityState[k]
					/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			starValues.c.v[0].val()[i] += starValues.c.v[k].cval()[i];
		}
		starValues.rho.val()[i] = densityState[0];
		starValues.v.val()[i] = velocityState;
		starValues.p.val()[i] = pressureState;
		starValues.a.val()[i] = aState;
		starValues.b.val()[i] = bState;
	}

	return std::make_tuple(numFluxes, starValues);
}
