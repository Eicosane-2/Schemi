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
	auto & mesh { surfaceOwnerSide.pressure.meshRef() };

	conservativeFlows numFluxes(mesh,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());
	starFields starValues(mesh,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());

	errors errValue { errors::noErrors };

	const std::valarray<scalar> sonicSpeedOwner(
			std::sqrt(
					surfaceOwnerSide.phaseThermodynamics->sqSonicSpeed(
							surfaceOwnerSide.concentration.p,
							surfaceOwnerSide.density[0].ref(),
							surfaceOwnerSide.internalEnergy.ref(),
							surfaceOwnerSide.pressure.ref())
							+ twothirds * surfaceOwnerSide.kTurb.ref()
									* (scalar(1.)
											+ surfaceOwnerSide.phaseThermodynamics->dpdUv(
													surfaceOwnerSide.concentration.p,
													surfaceOwnerSide.internalEnergy.ref()))));

	const std::valarray<scalar> sonicSpeedNeighbour(
			std::sqrt(
					surfaceNeighbourSide.phaseThermodynamics->sqSonicSpeed(
							surfaceNeighbourSide.concentration.p,
							surfaceNeighbourSide.density[0].ref(),
							surfaceNeighbourSide.internalEnergy.ref(),
							surfaceNeighbourSide.pressure.ref())
							+ twothirds * surfaceNeighbourSide.kTurb.ref()
									* (scalar(1.)
											+ surfaceNeighbourSide.phaseThermodynamics->dpdUv(
													surfaceNeighbourSide.concentration.p,
													surfaceNeighbourSide.internalEnergy.ref()))));

	for (std::size_t i = 0; i < mesh.surfacesSize(); ++i)
	{
		const scalar velocityProjectionOwner {
				surfaceOwnerSide.velocity.ref()[i] & mesh.surfaces()[i].N() };
		const scalar velocityProjectionNeighbour {
				surfaceNeighbourSide.velocity.ref()[i] & mesh.surfaces()[i].N() };

		const scalar turbulentPressureOwner { twothirds
				* surfaceOwnerSide.rhokTurb.ref()[i] };
		const scalar turbulentPressureNeighbour { twothirds
				* surfaceNeighbourSide.rhokTurb.ref()[i] };

		const scalar SOwner { std::min(
				velocityProjectionOwner - sonicSpeedOwner[i],
				velocityProjectionNeighbour - sonicSpeedNeighbour[i]) };
		const scalar SNeighbour { std::max(
				velocityProjectionOwner + sonicSpeedOwner[i],
				velocityProjectionNeighbour + sonicSpeedNeighbour[i]) };

		const scalar Sc { ((surfaceNeighbourSide.pressure.ref()[i]
				+ turbulentPressureNeighbour)
				- (surfaceOwnerSide.pressure.ref()[i] + turbulentPressureOwner)
				+ surfaceOwnerSide.density[0].ref()[i] * velocityProjectionOwner
						* (SOwner - velocityProjectionOwner)
				- surfaceNeighbourSide.density[0].ref()[i]
						* velocityProjectionNeighbour
						* (SNeighbour - velocityProjectionNeighbour))
				/ (surfaceOwnerSide.density[0].ref()[i]
						* (SOwner - velocityProjectionOwner)
						- surfaceNeighbourSide.density[0].ref()[i]
								* (SNeighbour - velocityProjectionNeighbour)) };

		const scalar POwnNei { 0.5
				* ((surfaceOwnerSide.pressure.ref()[i] + turbulentPressureOwner)
						+ (surfaceNeighbourSide.pressure.ref()[i]
								+ turbulentPressureNeighbour)
						+ surfaceOwnerSide.density[0].ref()[i]
								* (SOwner - velocityProjectionOwner)
								* (Sc - velocityProjectionOwner)
						+ surfaceNeighbourSide.density[0].ref()[i]
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
				densityState[k] = surfaceOwnerSide.density[k].ref()[i];

			momentumState = surfaceOwnerSide.momentum.ref()[i];
			totalEnergyState = surfaceOwnerSide.totalEnergy.ref()[i];
			rhokState = surfaceOwnerSide.rhokTurb.ref()[i];
			rhoEpsState = surfaceOwnerSide.rhoepsTurb.ref()[i];
			rhoaState = surfaceOwnerSide.rhoaTurb.ref()[i];
			rhobState = surfaceOwnerSide.rhobTurb.ref()[i];

			velocityState = surfaceOwnerSide.velocity.ref()[i];
			pressureState = surfaceOwnerSide.pressure.ref()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];
		}
		else if (Sc >= 0)
		{
			const scalar factor { SOwner - velocityProjectionOwner };
			const scalar denominator { 1 / (SOwner - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceOwnerSide.density[k].ref()[i] * factor
						* denominator;

			momentumState = (surfaceOwnerSide.momentum.ref()[i] * factor
					- mesh.surfaces()[i].N()
							* (surfaceOwnerSide.pressure.ref()[i]
									+ turbulentPressureOwner - POwnNei))
					* denominator;
			totalEnergyState = (surfaceOwnerSide.totalEnergy.ref()[i] * factor
					- (surfaceOwnerSide.pressure.ref()[i]
							+ turbulentPressureOwner) * velocityProjectionOwner
					+ POwnNei * Sc) * denominator;
			rhokState = surfaceOwnerSide.rhokTurb.ref()[i] * factor
					* denominator;
			rhoEpsState = surfaceOwnerSide.rhoepsTurb.ref()[i] * factor
					* denominator;
			rhoaState = surfaceOwnerSide.rhoaTurb.ref()[i] * factor
					* denominator;
			rhobState = surfaceOwnerSide.rhobTurb.ref()[i] * factor
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
				densityState[k] = surfaceNeighbourSide.density[k].ref()[i]
						* factor * denominator;

			momentumState = (surfaceNeighbourSide.momentum.ref()[i] * factor
					- mesh.surfaces()[i].N()
							* (surfaceNeighbourSide.pressure.ref()[i]
									+ turbulentPressureNeighbour - POwnNei))
					* denominator;
			totalEnergyState = (surfaceNeighbourSide.totalEnergy.ref()[i]
					* factor
					- (surfaceNeighbourSide.pressure.ref()[i]
							+ turbulentPressureNeighbour)
							* velocityProjectionNeighbour + POwnNei * Sc)
					* denominator;
			rhokState = surfaceNeighbourSide.rhokTurb.ref()[i] * factor
					* denominator;
			rhoEpsState = surfaceNeighbourSide.rhoepsTurb.ref()[i] * factor
					* denominator;
			rhoaState = surfaceNeighbourSide.rhoaTurb.ref()[i] * factor
					* denominator;
			rhobState = surfaceNeighbourSide.rhobTurb.ref()[i] * factor
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
				densityState[k] = surfaceNeighbourSide.density[k].ref()[i];

			momentumState = surfaceNeighbourSide.momentum.ref()[i];
			totalEnergyState = surfaceNeighbourSide.totalEnergy.ref()[i];
			rhokState = surfaceNeighbourSide.rhokTurb.ref()[i];
			rhoEpsState = surfaceNeighbourSide.rhoepsTurb.ref()[i];
			rhoaState = surfaceNeighbourSide.rhoaTurb.ref()[i];
			rhobState = surfaceNeighbourSide.rhobTurb.ref()[i];

			velocityState = surfaceNeighbourSide.velocity.ref()[i];
			pressureState = surfaceNeighbourSide.pressure.ref()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];
		}
		else
		{
			errValue = errors::RiemannSolverError;
			throw exception("Fatal error in Riemann solver.", errValue);
		}

		for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
			numFluxes.density[k].ref_r()[i] = velocityState * densityState[k];

		numFluxes.momentum.ref_r()[i] = momentumState * velocityState
				+ tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
						* (pressureState + twothirds * rhokState);

		numFluxes.totalEnergy.ref_r()[i] = velocityState * totalEnergyState
				+ velocityState * (pressureState + twothirds * rhokState);

		numFluxes.rhokTurb.ref_r()[i] = rhokState * velocityState;
		numFluxes.rhoepsTurb.ref_r()[i] = velocityState * rhoEpsState;
		numFluxes.rhoaTurb.ref_r()[i] = rhoaState * velocityState;
		numFluxes.rhobTurb.ref_r()[i] = velocityState * rhobState;

		starValues.c.v[0].ref_r()[i] = 0;
		for (std::size_t k = 1; k < densityState.size(); ++k)
		{
			starValues.c.v[k].ref_r()[i] = densityState[k]
					/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			starValues.c.v[0].ref_r()[i] += starValues.c.v[k].ref()[i];
		}
		starValues.rho.ref_r()[i] = densityState[0];
		starValues.v.ref_r()[i] = velocityState;
		starValues.p.ref_r()[i] = pressureState;
		starValues.a.ref_r()[i] = aState;
		starValues.b.ref_r()[i] = bState;
	}

	return std::make_tuple(numFluxes, starValues);
}
