/*
 * HLLC2pSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HLLC2pSolver.hpp"

#include "tensorVectorDotProduct.hpp"
#include "vectorVectorDotProduct.hpp"

std::tuple<schemi::conservativeFlows, schemi::starFields> schemi::HLLC2pSolver::calculateFlows(
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

		const scalar POwn { surfaceOwnerSide.pressure()[i]
				+ surfaceOwnerSide.density[0]()[i]
						* (SOwner - velocityProjectionOwner)
						* (Sc - velocityProjectionOwner) };

		const scalar PNei { surfaceNeighbourSide.pressure()[i]
				+ surfaceNeighbourSide.density[0]()[i]
						* (SNeighbour - velocityProjectionNeighbour)
						* (Sc - velocityProjectionNeighbour) };

		std::valarray<scalar> densityState(numFluxes.density.size());
		vector velocityState;
		scalar pressureState;
		vector aState;
		scalar bState;

		if (SOwner >= 0)
		{
			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceOwnerSide.density[k]()[i];

			const vector momentumState = surfaceOwnerSide.momentum()[i];
			const scalar totalEnergyState = surfaceOwnerSide.totalEnergy()[i];
			const scalar rhokState = surfaceOwnerSide.rhokTurb()[i];
			const scalar rhoEpsState = surfaceOwnerSide.rhoepsTurb()[i];
			const vector rhoaState = surfaceOwnerSide.rhoaTurb()[i];
			const scalar rhobState = surfaceOwnerSide.rhobTurb()[i];

			velocityState = surfaceOwnerSide.velocity()[i];
			pressureState = surfaceOwnerSide.pressure()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].r()[i] = velocityProjectionOwner
						* mesh_.surfaces()[i].N() * densityState[k];

			numFluxes.momentum.r()[i] = (momentumState * velocityProjectionOwner
					+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh_.surfaces()[i].N())
							* (pressureState + turbulentPressureOwner))
					* mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.r()[i] =
					(velocityProjectionOwner
							* (totalEnergyState + pressureState
									+ turbulentPressureOwner))
							* mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.r()[i] = velocityProjectionOwner * rhokState
					* mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.r()[i] = velocityProjectionOwner
					* mesh_.surfaces()[i].N() * rhoEpsState;
			numFluxes.rhoaTurb.r()[i] = rhoaState * velocityProjectionOwner
					* mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.r()[i] = velocityProjectionOwner
					* mesh_.surfaces()[i].N() * rhobState;
		}
		else if (Sc >= 0)
		{
			const scalar factor { SOwner - velocityProjectionOwner };
			const scalar denominator { 1 / (SOwner - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceOwnerSide.density[k]()[i] * factor
						* denominator;

			const vector momentumState = (surfaceOwnerSide.momentum()[i]
					* factor
					- mesh_.surfaces()[i].N()
							* (surfaceOwnerSide.pressure()[i]
									+ turbulentPressureOwner - POwn))
					* denominator;
			const scalar totalEnergyState = (surfaceOwnerSide.totalEnergy()[i]
					* factor
					- (surfaceOwnerSide.pressure()[i] + turbulentPressureOwner)
							* velocityProjectionOwner + POwn * Sc)
					* denominator;
			const scalar rhokState = surfaceOwnerSide.rhokTurb()[i] * factor
					* denominator;
			const scalar rhoEpsState = surfaceOwnerSide.rhoepsTurb()[i] * factor
					* denominator;
			const vector rhoaState = surfaceOwnerSide.rhoaTurb()[i] * factor
					* denominator;
			const scalar rhobState = surfaceOwnerSide.rhobTurb()[i] * factor
					* denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceNeighbourSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].r()[i] = (surfaceOwnerSide.density[k]()[i]
						* velocityProjectionOwner
						+ SOwner
								* (densityState[k]
										- surfaceOwnerSide.density[k]()[i]))
						* mesh_.surfaces()[i].N();

			numFluxes.momentum.r()[i] = (surfaceOwnerSide.momentum()[i]
					* velocityProjectionOwner
					+ ((tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh_.surfaces()[i].N())
							* (surfaceOwnerSide.pressure()[i]
									+ turbulentPressureOwner))
					+ SOwner * (momentumState - surfaceOwnerSide.momentum()[i]))
					* mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.r()[i] = ((surfaceOwnerSide.totalEnergy()[i]
					+ surfaceOwnerSide.pressure()[i] + turbulentPressureOwner)
					* velocityProjectionOwner
					+ SOwner
							* (totalEnergyState
									- surfaceOwnerSide.totalEnergy()[i]))
					* mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.r()[i] = (surfaceOwnerSide.rhokTurb()[i]
					* velocityProjectionOwner
					+ SOwner * (rhokState - surfaceOwnerSide.rhokTurb()[i]))
					* mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.r()[i] = (surfaceOwnerSide.rhoepsTurb()[i]
					* velocityProjectionOwner
					+ SOwner * (rhoEpsState - surfaceOwnerSide.rhoepsTurb()[i]))
					* mesh_.surfaces()[i].N();
			numFluxes.rhoaTurb.r()[i] = (surfaceOwnerSide.rhoaTurb()[i]
					* velocityProjectionOwner
					+ SOwner * (rhoaState - surfaceOwnerSide.rhoaTurb()[i]))
					* mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.r()[i] = (surfaceOwnerSide.rhobTurb()[i]
					* velocityProjectionOwner
					+ SOwner * (rhobState - surfaceOwnerSide.rhobTurb()[i]))
					* mesh_.surfaces()[i].N();
		}
		else if (SNeighbour > 0)
		{
			const scalar factor { SNeighbour - velocityProjectionNeighbour };
			const scalar denominator { 1 / (SNeighbour - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceNeighbourSide.density[k]()[i] * factor
						* denominator;

			const vector momentumState = (surfaceNeighbourSide.momentum()[i]
					* factor
					- mesh_.surfaces()[i].N()
							* (surfaceNeighbourSide.pressure()[i]
									+ turbulentPressureNeighbour - PNei))
					* denominator;
			const scalar totalEnergyState =
					(surfaceNeighbourSide.totalEnergy()[i] * factor
							- (surfaceNeighbourSide.pressure()[i]
									+ turbulentPressureNeighbour)
									* velocityProjectionNeighbour + PNei * Sc)
							* denominator;
			const scalar rhokState = surfaceNeighbourSide.rhokTurb()[i] * factor
					* denominator;
			const scalar rhoEpsState = surfaceNeighbourSide.rhoepsTurb()[i]
					* factor * denominator;
			const vector rhoaState = surfaceNeighbourSide.rhoaTurb()[i] * factor
					* denominator;
			const scalar rhobState = surfaceNeighbourSide.rhobTurb()[i] * factor
					* denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceNeighbourSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].r()[i] =
						(surfaceNeighbourSide.density[k]()[i]
								* velocityProjectionNeighbour
								+ SNeighbour
										* (densityState[k]
												- surfaceNeighbourSide.density[k]()[i]))
								* mesh_.surfaces()[i].N();

			numFluxes.momentum.r()[i] = (surfaceNeighbourSide.momentum()[i]
					* velocityProjectionNeighbour
					+ ((tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh_.surfaces()[i].N())
							* (surfaceNeighbourSide.pressure()[i]
									+ turbulentPressureNeighbour))
					+ SNeighbour
							* (momentumState
									- surfaceNeighbourSide.momentum()[i]))
					* mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.r()[i] =
					((surfaceNeighbourSide.totalEnergy()[i]
							+ surfaceNeighbourSide.pressure()[i]
							+ turbulentPressureNeighbour)
							* velocityProjectionNeighbour
							+ SNeighbour
									* (totalEnergyState
											- surfaceNeighbourSide.totalEnergy()[i]))
							* mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.r()[i] = (surfaceNeighbourSide.rhokTurb()[i]
					* velocityProjectionNeighbour
					+ SNeighbour
							* (rhokState - surfaceNeighbourSide.rhokTurb()[i]))
					* mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.r()[i] = (surfaceNeighbourSide.rhoepsTurb()[i]
					* velocityProjectionNeighbour
					+ SNeighbour
							* (rhoEpsState
									- surfaceNeighbourSide.rhoepsTurb()[i]))
					* mesh_.surfaces()[i].N();
			numFluxes.rhoaTurb.r()[i] = (surfaceNeighbourSide.rhoaTurb()[i]
					* velocityProjectionNeighbour
					+ SNeighbour
							* (rhoaState - surfaceNeighbourSide.rhoaTurb()[i]))
					* mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.r()[i] = (surfaceNeighbourSide.rhobTurb()[i]
					* velocityProjectionNeighbour
					+ SNeighbour
							* (rhobState - surfaceNeighbourSide.rhobTurb()[i]))
					* mesh_.surfaces()[i].N();
		}
		else if (SNeighbour <= 0)
		{
			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceNeighbourSide.density[k]()[i];

			const vector momentumState = surfaceNeighbourSide.momentum()[i];
			const scalar totalEnergyState =
					surfaceNeighbourSide.totalEnergy()[i];
			const scalar rhokState = surfaceNeighbourSide.rhokTurb()[i];
			const scalar rhoEpsState = surfaceNeighbourSide.rhoepsTurb()[i];
			const vector rhoaState = surfaceNeighbourSide.rhoaTurb()[i];
			const scalar rhobState = surfaceNeighbourSide.rhobTurb()[i];

			velocityState = surfaceNeighbourSide.velocity()[i];
			pressureState = surfaceNeighbourSide.pressure()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].r()[i] = velocityProjectionNeighbour
						* mesh_.surfaces()[i].N() * densityState[k];

			numFluxes.momentum.r()[i] = (momentumState
					* velocityProjectionNeighbour
					+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh_.surfaces()[i].N())
							* (pressureState + turbulentPressureNeighbour))
					* mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.r()[i] = velocityProjectionNeighbour
					* (totalEnergyState + pressureState
							+ turbulentPressureNeighbour)
					* mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.r()[i] = velocityProjectionNeighbour * rhokState
					* mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.r()[i] = velocityProjectionNeighbour
					* mesh_.surfaces()[i].N() * rhoEpsState;
			numFluxes.rhoaTurb.r()[i] = rhoaState * velocityProjectionNeighbour
					* mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.r()[i] = velocityProjectionNeighbour
					* mesh_.surfaces()[i].N() * rhobState;
		}
		else
		{
			errValue = errors::RiemannSolverError;
			throw exception("Fatal error in Riemann solver.", errValue);
		}

		starValues.c.v[0].r()[i] = 0;
		for (std::size_t k = 1; k < densityState.size(); ++k)
		{
			starValues.c.v[k].r()[i] = densityState[k]
					/ surfaceNeighbourSide.phaseThermodynamics->Mv()[k - 1];

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

