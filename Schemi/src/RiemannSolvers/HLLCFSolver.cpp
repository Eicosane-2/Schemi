/*
 * HLLCFSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HLLCFSolver.hpp"

#include "vectorVectorDotProduct.hpp"

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
							surfaceOwnerSide.density[0](),
							surfaceOwnerSide.internalEnergy(),
							surfaceOwnerSide.pressure())
							+ twothirds * surfaceOwnerSide.kTurb()
									* (scalar(1.)
											+ surfaceOwnerSide.phaseThermodynamics->dpdUv(
													surfaceOwnerSide.concentration.p,
													surfaceOwnerSide.internalEnergy()))));

	const std::valarray<scalar> sonicSpeedNeighbour(
			std::sqrt(
					surfaceNeighbourSide.phaseThermodynamics->sqSonicSpeed(
							surfaceNeighbourSide.concentration.p,
							surfaceNeighbourSide.density[0](),
							surfaceNeighbourSide.internalEnergy(),
							surfaceNeighbourSide.pressure())
							+ twothirds * surfaceNeighbourSide.kTurb()
									* (scalar(1.)
											+ surfaceNeighbourSide.phaseThermodynamics->dpdUv(
													surfaceNeighbourSide.concentration.p,
													surfaceNeighbourSide.internalEnergy()))));

	for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
	{
		const scalar velocityProjectionOwner { surfaceOwnerSide.velocity()[i]
				& mesh_.surfaces()[i].N() };
		const scalar velocityProjectionNeighbour {
				surfaceNeighbourSide.velocity()[i] & mesh_.surfaces()[i].N() };

		const scalar turbulentPressureOwner { twothirds
				* surfaceOwnerSide.rhokTurb()[i] };
		const scalar turbulentPressureNeighbour { twothirds
				* surfaceNeighbourSide.rhokTurb()[i] };

		const scalar SOwner { std::min(
				velocityProjectionOwner - sonicSpeedOwner[i],
				velocityProjectionNeighbour - sonicSpeedNeighbour[i]) };
		const scalar SNeighbour { std::max(
				velocityProjectionOwner + sonicSpeedOwner[i],
				velocityProjectionNeighbour + sonicSpeedNeighbour[i]) };

		const scalar Sc { ((surfaceNeighbourSide.pressure()[i]
				+ turbulentPressureNeighbour)
				- (surfaceOwnerSide.pressure()[i] + turbulentPressureOwner)
				+ surfaceOwnerSide.density[0]()[i] * velocityProjectionOwner
						* (SOwner - velocityProjectionOwner)
				- surfaceNeighbourSide.density[0]()[i]
						* velocityProjectionNeighbour
						* (SNeighbour - velocityProjectionNeighbour))
				/ (surfaceOwnerSide.density[0]()[i]
						* (SOwner - velocityProjectionOwner)
						- surfaceNeighbourSide.density[0]()[i]
								* (SNeighbour - velocityProjectionNeighbour)) };

		const scalar POwnNei { 0.5
				* ((surfaceOwnerSide.pressure()[i] + turbulentPressureOwner)
						+ (surfaceNeighbourSide.pressure()[i]
								+ turbulentPressureNeighbour)
						+ surfaceOwnerSide.density[0]()[i]
								* (SOwner - velocityProjectionOwner)
								* (Sc - velocityProjectionOwner)
						+ surfaceNeighbourSide.density[0]()[i]
								* (SNeighbour - velocityProjectionNeighbour)
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
				densityState[k] = surfaceOwnerSide.density[k]()[i];

			momentumState = surfaceOwnerSide.momentum()[i];
			totalEnergyState = surfaceOwnerSide.totalEnergy()[i];
			rhokState = surfaceOwnerSide.rhokTurb()[i];
			rhoEpsState = surfaceOwnerSide.rhoepsTurb()[i];
			rhoaState = surfaceOwnerSide.rhoaTurb()[i];
			rhobState = surfaceOwnerSide.rhobTurb()[i];

			velocityState = surfaceOwnerSide.velocity()[i];
			pressureState = surfaceOwnerSide.pressure()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];
		}
		else if (Sc >= 0)
		{
			const scalar factor { SOwner - velocityProjectionOwner };
			const scalar denominator { 1 / (SOwner - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceOwnerSide.density[k]()[i] * factor
						* denominator;

			momentumState = (surfaceOwnerSide.momentum()[i] * factor
					- mesh_.surfaces()[i].N()
							* (surfaceOwnerSide.pressure()[i]
									+ turbulentPressureOwner - POwnNei))
					* denominator;
			totalEnergyState = (surfaceOwnerSide.totalEnergy()[i] * factor
					- (surfaceOwnerSide.pressure()[i] + turbulentPressureOwner)
							* velocityProjectionOwner + POwnNei * Sc)
					* denominator;
			rhokState = surfaceOwnerSide.rhokTurb()[i] * factor * denominator;
			rhoEpsState = surfaceOwnerSide.rhoepsTurb()[i] * factor
					* denominator;
			rhoaState = surfaceOwnerSide.rhoaTurb()[i] * factor * denominator;
			rhobState = surfaceOwnerSide.rhobTurb()[i] * factor * denominator;

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
				densityState[k] = surfaceNeighbourSide.density[k]()[i] * factor
						* denominator;

			momentumState = (surfaceNeighbourSide.momentum()[i] * factor
					- mesh_.surfaces()[i].N()
							* (surfaceNeighbourSide.pressure()[i]
									+ turbulentPressureNeighbour - POwnNei))
					* denominator;
			totalEnergyState = (surfaceNeighbourSide.totalEnergy()[i] * factor
					- (surfaceNeighbourSide.pressure()[i]
							+ turbulentPressureNeighbour)
							* velocityProjectionNeighbour + POwnNei * Sc)
					* denominator;
			rhokState = surfaceNeighbourSide.rhokTurb()[i] * factor
					* denominator;
			rhoEpsState = surfaceNeighbourSide.rhoepsTurb()[i] * factor
					* denominator;
			rhoaState = surfaceNeighbourSide.rhoaTurb()[i] * factor
					* denominator;
			rhobState = surfaceNeighbourSide.rhobTurb()[i] * factor
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
				densityState[k] = surfaceNeighbourSide.density[k]()[i];

			momentumState = surfaceNeighbourSide.momentum()[i];
			totalEnergyState = surfaceNeighbourSide.totalEnergy()[i];
			rhokState = surfaceNeighbourSide.rhokTurb()[i];
			rhoEpsState = surfaceNeighbourSide.rhoepsTurb()[i];
			rhoaState = surfaceNeighbourSide.rhoaTurb()[i];
			rhobState = surfaceNeighbourSide.rhobTurb()[i];

			velocityState = surfaceNeighbourSide.velocity()[i];
			pressureState = surfaceNeighbourSide.pressure()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];
		}
		else
		{
			errValue = errors::RiemannSolverError;
			throw exception("Fatal error in Riemann solver.", errValue);
		}

		for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
			numFluxes.density[k].r()[i] = velocityState * densityState[k];

		numFluxes.momentum.r()[i] = momentumState * velocityState
				+ tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
						* (pressureState + twothirds * rhokState);

		numFluxes.totalEnergy.r()[i] = velocityState * totalEnergyState
				+ velocityState * (pressureState + twothirds * rhokState);

		numFluxes.rhokTurb.r()[i] = rhokState * velocityState;
		numFluxes.rhoepsTurb.r()[i] = velocityState * rhoEpsState;
		numFluxes.rhoaTurb.r()[i] = rhoaState * velocityState;
		numFluxes.rhobTurb.r()[i] = velocityState * rhobState;

		starValues.c.v[0].r()[i] = 0;
		for (std::size_t k = 1; k < densityState.size(); ++k)
		{
			starValues.c.v[k].r()[i] = densityState[k]
					/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			starValues.c.v[0].r()[i] += starValues.c.v[k]()[i];
		}
		starValues.rho.r()[i] = densityState[0];
		starValues.v.r()[i] = velocityState;
		starValues.p.r()[i] = pressureState;
		starValues.a.r()[i] = aState;
		starValues.b.r()[i] = bState;
	}

	return std::make_tuple(numFluxes, starValues);
}
