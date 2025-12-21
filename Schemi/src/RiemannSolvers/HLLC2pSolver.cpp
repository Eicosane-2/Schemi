/*
 * HLLC2pSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HLLC2pSolver.hpp"

#include "vector.hpp"
#include "tensor.hpp"

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

		const scalar POwn { surfaceOwnerSide.pressure.cval()[i]
				+ surfaceOwnerSide.density[0].cval()[i]
						* (SOwner - velocityProjectionOwner)
						* (Sc - velocityProjectionOwner) };

		const scalar PNei { surfaceNeighbourSide.pressure.cval()[i]
				+ surfaceNeighbourSide.density[0].cval()[i]
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
				densityState[k] = surfaceOwnerSide.density[k].cval()[i];

			const vector momentumState = surfaceOwnerSide.momentum.cval()[i];
			const scalar totalEnergyState =
					surfaceOwnerSide.totalEnergy.cval()[i];
			const scalar rhokState = surfaceOwnerSide.rhokTurb.cval()[i];
			const scalar rhoEpsState = surfaceOwnerSide.rhoepsTurb.cval()[i];
			const vector rhoaState = surfaceOwnerSide.rhoaTurb.cval()[i];
			const scalar rhobState = surfaceOwnerSide.rhobTurb.cval()[i];

			velocityState = surfaceOwnerSide.velocity.cval()[i];
			pressureState = surfaceOwnerSide.pressure.cval()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].val()[i] = velocityProjectionOwner
						* mesh_.surfaces()[i].N() * densityState[k];

			numFluxes.momentum.val()[i] = (momentumState
					* velocityProjectionOwner
					+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh_.surfaces()[i].N())
							* (pressureState + turbulentPressureOwner))
					* mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.val()[i] =
					(velocityProjectionOwner
							* (totalEnergyState + pressureState
									+ turbulentPressureOwner))
							* mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.val()[i] = velocityProjectionOwner * rhokState
					* mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.val()[i] = velocityProjectionOwner
					* mesh_.surfaces()[i].N() * rhoEpsState;
			numFluxes.rhoaTurb.val()[i] = rhoaState * velocityProjectionOwner
					* mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.val()[i] = velocityProjectionOwner
					* mesh_.surfaces()[i].N() * rhobState;
		}
		else if (Sc >= 0)
		{
			const scalar factor { SOwner - velocityProjectionOwner };
			const scalar denominator { 1 / (SOwner - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceOwnerSide.density[k].cval()[i] * factor
						* denominator;

			const vector momentumState = (surfaceOwnerSide.momentum.cval()[i]
					* factor
					- mesh_.surfaces()[i].N()
							* (surfaceOwnerSide.pressure.cval()[i]
									+ turbulentPressureOwner - POwn))
					* denominator;
			const scalar totalEnergyState =
					(surfaceOwnerSide.totalEnergy.cval()[i] * factor
							- (surfaceOwnerSide.pressure.cval()[i]
									+ turbulentPressureOwner)
									* velocityProjectionOwner + POwn * Sc)
							* denominator;
			const scalar rhokState = surfaceOwnerSide.rhokTurb.cval()[i]
					* factor * denominator;
			const scalar rhoEpsState = surfaceOwnerSide.rhoepsTurb.cval()[i]
					* factor * denominator;
			const vector rhoaState = surfaceOwnerSide.rhoaTurb.cval()[i]
					* factor * denominator;
			const scalar rhobState = surfaceOwnerSide.rhobTurb.cval()[i]
					* factor * denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceNeighbourSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].val()[i] =
						(surfaceOwnerSide.density[k].cval()[i]
								* velocityProjectionOwner
								+ SOwner
										* (densityState[k]
												- surfaceOwnerSide.density[k].cval()[i]))
								* mesh_.surfaces()[i].N();

			numFluxes.momentum.val()[i] = (surfaceOwnerSide.momentum.cval()[i]
					* velocityProjectionOwner
					+ ((tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh_.surfaces()[i].N())
							* (surfaceOwnerSide.pressure.cval()[i]
									+ turbulentPressureOwner))
					+ SOwner
							* (momentumState
									- surfaceOwnerSide.momentum.cval()[i]))
					* mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.val()[i] =
					((surfaceOwnerSide.totalEnergy.cval()[i]
							+ surfaceOwnerSide.pressure.cval()[i]
							+ turbulentPressureOwner) * velocityProjectionOwner
							+ SOwner
									* (totalEnergyState
											- surfaceOwnerSide.totalEnergy.cval()[i]))
							* mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.val()[i] =
					(surfaceOwnerSide.rhokTurb.cval()[i]
							* velocityProjectionOwner
							+ SOwner
									* (rhokState
											- surfaceOwnerSide.rhokTurb.cval()[i]))
							* mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.val()[i] =
					(surfaceOwnerSide.rhoepsTurb.cval()[i]
							* velocityProjectionOwner
							+ SOwner
									* (rhoEpsState
											- surfaceOwnerSide.rhoepsTurb.cval()[i]))
							* mesh_.surfaces()[i].N();
			numFluxes.rhoaTurb.val()[i] =
					(surfaceOwnerSide.rhoaTurb.cval()[i]
							* velocityProjectionOwner
							+ SOwner
									* (rhoaState
											- surfaceOwnerSide.rhoaTurb.cval()[i]))
							* mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.val()[i] =
					(surfaceOwnerSide.rhobTurb.cval()[i]
							* velocityProjectionOwner
							+ SOwner
									* (rhobState
											- surfaceOwnerSide.rhobTurb.cval()[i]))
							* mesh_.surfaces()[i].N();
		}
		else if (SNeighbour > 0)
		{
			const scalar factor { SNeighbour - velocityProjectionNeighbour };
			const scalar denominator { 1 / (SNeighbour - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceNeighbourSide.density[k].cval()[i]
						* factor * denominator;

			const vector momentumState =
					(surfaceNeighbourSide.momentum.cval()[i] * factor
							- mesh_.surfaces()[i].N()
									* (surfaceNeighbourSide.pressure.cval()[i]
											+ turbulentPressureNeighbour - PNei))
							* denominator;
			const scalar totalEnergyState =
					(surfaceNeighbourSide.totalEnergy.cval()[i] * factor
							- (surfaceNeighbourSide.pressure.cval()[i]
									+ turbulentPressureNeighbour)
									* velocityProjectionNeighbour + PNei * Sc)
							* denominator;
			const scalar rhokState = surfaceNeighbourSide.rhokTurb.cval()[i]
					* factor * denominator;
			const scalar rhoEpsState = surfaceNeighbourSide.rhoepsTurb.cval()[i]
					* factor * denominator;
			const vector rhoaState = surfaceNeighbourSide.rhoaTurb.cval()[i]
					* factor * denominator;
			const scalar rhobState = surfaceNeighbourSide.rhobTurb.cval()[i]
					* factor * denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceNeighbourSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].val()[i] =
						(surfaceNeighbourSide.density[k].cval()[i]
								* velocityProjectionNeighbour
								+ SNeighbour
										* (densityState[k]
												- surfaceNeighbourSide.density[k].cval()[i]))
								* mesh_.surfaces()[i].N();

			numFluxes.momentum.val()[i] =
					(surfaceNeighbourSide.momentum.cval()[i]
							* velocityProjectionNeighbour
							+ ((tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
									& mesh_.surfaces()[i].N())
									* (surfaceNeighbourSide.pressure.cval()[i]
											+ turbulentPressureNeighbour))
							+ SNeighbour
									* (momentumState
											- surfaceNeighbourSide.momentum.cval()[i]))
							* mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.val()[i] =
					((surfaceNeighbourSide.totalEnergy.cval()[i]
							+ surfaceNeighbourSide.pressure.cval()[i]
							+ turbulentPressureNeighbour)
							* velocityProjectionNeighbour
							+ SNeighbour
									* (totalEnergyState
											- surfaceNeighbourSide.totalEnergy.cval()[i]))
							* mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.val()[i] =
					(surfaceNeighbourSide.rhokTurb.cval()[i]
							* velocityProjectionNeighbour
							+ SNeighbour
									* (rhokState
											- surfaceNeighbourSide.rhokTurb.cval()[i]))
							* mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.val()[i] =
					(surfaceNeighbourSide.rhoepsTurb.cval()[i]
							* velocityProjectionNeighbour
							+ SNeighbour
									* (rhoEpsState
											- surfaceNeighbourSide.rhoepsTurb.cval()[i]))
							* mesh_.surfaces()[i].N();
			numFluxes.rhoaTurb.val()[i] =
					(surfaceNeighbourSide.rhoaTurb.cval()[i]
							* velocityProjectionNeighbour
							+ SNeighbour
									* (rhoaState
											- surfaceNeighbourSide.rhoaTurb.cval()[i]))
							* mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.val()[i] =
					(surfaceNeighbourSide.rhobTurb.cval()[i]
							* velocityProjectionNeighbour
							+ SNeighbour
									* (rhobState
											- surfaceNeighbourSide.rhobTurb.cval()[i]))
							* mesh_.surfaces()[i].N();
		}
		else if (SNeighbour <= 0)
		{
			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceNeighbourSide.density[k].cval()[i];

			const vector momentumState = surfaceNeighbourSide.momentum.cval()[i];
			const scalar totalEnergyState =
					surfaceNeighbourSide.totalEnergy.cval()[i];
			const scalar rhokState = surfaceNeighbourSide.rhokTurb.cval()[i];
			const scalar rhoEpsState = surfaceNeighbourSide.rhoepsTurb.cval()[i];
			const vector rhoaState = surfaceNeighbourSide.rhoaTurb.cval()[i];
			const scalar rhobState = surfaceNeighbourSide.rhobTurb.cval()[i];

			velocityState = surfaceNeighbourSide.velocity.cval()[i];
			pressureState = surfaceNeighbourSide.pressure.cval()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].val()[i] = velocityProjectionNeighbour
						* mesh_.surfaces()[i].N() * densityState[k];

			numFluxes.momentum.val()[i] = (momentumState
					* velocityProjectionNeighbour
					+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh_.surfaces()[i].N())
							* (pressureState + turbulentPressureNeighbour))
					* mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.val()[i] = velocityProjectionNeighbour
					* (totalEnergyState + pressureState
							+ turbulentPressureNeighbour)
					* mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.val()[i] = velocityProjectionNeighbour
					* rhokState * mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.val()[i] = velocityProjectionNeighbour
					* mesh_.surfaces()[i].N() * rhoEpsState;
			numFluxes.rhoaTurb.val()[i] = rhoaState
					* velocityProjectionNeighbour * mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.val()[i] = velocityProjectionNeighbour
					* mesh_.surfaces()[i].N() * rhobState;
		}
		else
			[[unlikely]]
			throw exception("Fatal error in Riemann solver.", errValue =
					errors::RiemannSolverError);

		starValues.c.v[0].val()[i] = 0;
		for (std::size_t k = 1; k < densityState.size(); ++k)
		{
			starValues.c.v[k].val()[i] = densityState[k]
					/ surfaceNeighbourSide.phaseThermodynamics->Mv()[k - 1];

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

