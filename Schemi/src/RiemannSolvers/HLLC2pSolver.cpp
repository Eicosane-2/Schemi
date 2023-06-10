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
	auto & mesh { surfaceOwnerSide.pressure.meshRef() };

	conservativeFlows numFluxes(mesh,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());
	starFields starValues(mesh,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());

	errorsEnum errValue { errorsEnum::noErrors };

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

		const scalar POwn { surfaceOwnerSide.pressure.ref()[i]
				+ surfaceOwnerSide.density[0].ref()[i]
						* (SOwner - velocityProjectionOwner)
						* (Sc - velocityProjectionOwner) };

		const scalar PNei { surfaceNeighbourSide.pressure.ref()[i]
				+ surfaceNeighbourSide.density[0].ref()[i]
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
				densityState[k] = surfaceOwnerSide.density[k].ref()[i];

			const vector momentumState = surfaceOwnerSide.momentum.ref()[i];
			const scalar totalEnergyState =
					surfaceOwnerSide.totalEnergy.ref()[i];
			const scalar rhokState = surfaceOwnerSide.rhokTurb.ref()[i];
			const scalar rhoEpsState = surfaceOwnerSide.rhoepsTurb.ref()[i];
			const vector rhoaState = surfaceOwnerSide.rhoaTurb.ref()[i];
			const scalar rhobState = surfaceOwnerSide.rhobTurb.ref()[i];

			velocityState = surfaceOwnerSide.velocity.ref()[i];
			pressureState = surfaceOwnerSide.pressure.ref()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].ref_r()[i] = velocityProjectionOwner
						* mesh.surfaces()[i].N() * densityState[k];

			numFluxes.momentum.ref_r()[i] = (momentumState
					* velocityProjectionOwner
					+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh.surfaces()[i].N())
							* (pressureState + turbulentPressureOwner))
					* mesh.surfaces()[i].N();

			numFluxes.totalEnergy.ref_r()[i] =
					(velocityProjectionOwner
							* (totalEnergyState + pressureState
									+ turbulentPressureOwner))
							* mesh.surfaces()[i].N();

			numFluxes.rhokTurb.ref_r()[i] = velocityProjectionOwner * rhokState
					* mesh.surfaces()[i].N();
			numFluxes.rhoepsTurb.ref_r()[i] = velocityProjectionOwner
					* mesh.surfaces()[i].N() * rhoEpsState;
			numFluxes.rhoaTurb.ref_r()[i] = rhoaState * velocityProjectionOwner
					* mesh.surfaces()[i].N();
			numFluxes.rhobTurb.ref_r()[i] = velocityProjectionOwner
					* mesh.surfaces()[i].N() * rhobState;
		}
		else if (Sc >= 0)
		{
			const scalar factor { SOwner - velocityProjectionOwner };
			const scalar denominator { 1 / (SOwner - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceOwnerSide.density[k].ref()[i] * factor
						* denominator;

			const vector momentumState = (surfaceOwnerSide.momentum.ref()[i]
					* factor
					- mesh.surfaces()[i].N()
							* (surfaceOwnerSide.pressure.ref()[i]
									+ turbulentPressureOwner - POwn))
					* denominator;
			const scalar totalEnergyState =
					(surfaceOwnerSide.totalEnergy.ref()[i] * factor
							- (surfaceOwnerSide.pressure.ref()[i]
									+ turbulentPressureOwner)
									* velocityProjectionOwner + POwn * Sc)
							* denominator;
			const scalar rhokState = surfaceOwnerSide.rhokTurb.ref()[i] * factor
					* denominator;
			const scalar rhoEpsState = surfaceOwnerSide.rhoepsTurb.ref()[i]
					* factor * denominator;
			const vector rhoaState = surfaceOwnerSide.rhoaTurb.ref()[i] * factor
					* denominator;
			const scalar rhobState = surfaceOwnerSide.rhobTurb.ref()[i] * factor
					* denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceNeighbourSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].ref_r()[i] =
						(surfaceOwnerSide.density[k].ref()[i]
								* velocityProjectionOwner
								+ SOwner
										* (densityState[k]
												- surfaceOwnerSide.density[k].ref()[i]))
								* mesh.surfaces()[i].N();

			numFluxes.momentum.ref_r()[i] = (surfaceOwnerSide.momentum.ref()[i]
					* velocityProjectionOwner
					+ ((tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh.surfaces()[i].N())
							* (surfaceOwnerSide.pressure.ref()[i]
									+ turbulentPressureOwner))
					+ SOwner
							* (momentumState
									- surfaceOwnerSide.momentum.ref()[i]))
					* mesh.surfaces()[i].N();

			numFluxes.totalEnergy.ref_r()[i] =
					((surfaceOwnerSide.totalEnergy.ref()[i]
							+ surfaceOwnerSide.pressure.ref()[i]
							+ turbulentPressureOwner) * velocityProjectionOwner
							+ SOwner
									* (totalEnergyState
											- surfaceOwnerSide.totalEnergy.ref()[i]))
							* mesh.surfaces()[i].N();

			numFluxes.rhokTurb.ref_r()[i] = (surfaceOwnerSide.rhokTurb.ref()[i]
					* velocityProjectionOwner
					+ SOwner * (rhokState - surfaceOwnerSide.rhokTurb.ref()[i]))
					* mesh.surfaces()[i].N();
			numFluxes.rhoepsTurb.ref_r()[i] =
					(surfaceOwnerSide.rhoepsTurb.ref()[i]
							* velocityProjectionOwner
							+ SOwner
									* (rhoEpsState
											- surfaceOwnerSide.rhoepsTurb.ref()[i]))
							* mesh.surfaces()[i].N();
			numFluxes.rhoaTurb.ref_r()[i] = (surfaceOwnerSide.rhoaTurb.ref()[i]
					* velocityProjectionOwner
					+ SOwner * (rhoaState - surfaceOwnerSide.rhoaTurb.ref()[i]))
					* mesh.surfaces()[i].N();
			numFluxes.rhobTurb.ref_r()[i] = (surfaceOwnerSide.rhobTurb.ref()[i]
					* velocityProjectionOwner
					+ SOwner * (rhobState - surfaceOwnerSide.rhobTurb.ref()[i]))
					* mesh.surfaces()[i].N();
		}
		else if (SNeighbour > 0)
		{
			const scalar factor { SNeighbour - velocityProjectionNeighbour };
			const scalar denominator { 1 / (SNeighbour - Sc) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceNeighbourSide.density[k].ref()[i]
						* factor * denominator;

			const vector momentumState = (surfaceNeighbourSide.momentum.ref()[i]
					* factor
					- mesh.surfaces()[i].N()
							* (surfaceNeighbourSide.pressure.ref()[i]
									+ turbulentPressureNeighbour - PNei))
					* denominator;
			const scalar totalEnergyState =
					(surfaceNeighbourSide.totalEnergy.ref()[i] * factor
							- (surfaceNeighbourSide.pressure.ref()[i]
									+ turbulentPressureNeighbour)
									* velocityProjectionNeighbour + PNei * Sc)
							* denominator;
			const scalar rhokState = surfaceNeighbourSide.rhokTurb.ref()[i]
					* factor * denominator;
			const scalar rhoEpsState = surfaceNeighbourSide.rhoepsTurb.ref()[i]
					* factor * denominator;
			const vector rhoaState = surfaceNeighbourSide.rhoaTurb.ref()[i]
					* factor * denominator;
			const scalar rhobState = surfaceNeighbourSide.rhobTurb.ref()[i]
					* factor * denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceNeighbourSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].ref_r()[i] =
						(surfaceNeighbourSide.density[k].ref()[i]
								* velocityProjectionNeighbour
								+ SNeighbour
										* (densityState[k]
												- surfaceNeighbourSide.density[k].ref()[i]))
								* mesh.surfaces()[i].N();

			numFluxes.momentum.ref_r()[i] =
					(surfaceNeighbourSide.momentum.ref()[i]
							* velocityProjectionNeighbour
							+ ((tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
									& mesh.surfaces()[i].N())
									* (surfaceNeighbourSide.pressure.ref()[i]
											+ turbulentPressureNeighbour))
							+ SNeighbour
									* (momentumState
											- surfaceNeighbourSide.momentum.ref()[i]))
							* mesh.surfaces()[i].N();

			numFluxes.totalEnergy.ref_r()[i] =
					((surfaceNeighbourSide.totalEnergy.ref()[i]
							+ surfaceNeighbourSide.pressure.ref()[i]
							+ turbulentPressureNeighbour)
							* velocityProjectionNeighbour
							+ SNeighbour
									* (totalEnergyState
											- surfaceNeighbourSide.totalEnergy.ref()[i]))
							* mesh.surfaces()[i].N();

			numFluxes.rhokTurb.ref_r()[i] =
					(surfaceNeighbourSide.rhokTurb.ref()[i]
							* velocityProjectionNeighbour
							+ SNeighbour
									* (rhokState
											- surfaceNeighbourSide.rhokTurb.ref()[i]))
							* mesh.surfaces()[i].N();
			numFluxes.rhoepsTurb.ref_r()[i] =
					(surfaceNeighbourSide.rhoepsTurb.ref()[i]
							* velocityProjectionNeighbour
							+ SNeighbour
									* (rhoEpsState
											- surfaceNeighbourSide.rhoepsTurb.ref()[i]))
							* mesh.surfaces()[i].N();
			numFluxes.rhoaTurb.ref_r()[i] =
					(surfaceNeighbourSide.rhoaTurb.ref()[i]
							* velocityProjectionNeighbour
							+ SNeighbour
									* (rhoaState
											- surfaceNeighbourSide.rhoaTurb.ref()[i]))
							* mesh.surfaces()[i].N();
			numFluxes.rhobTurb.ref_r()[i] =
					(surfaceNeighbourSide.rhobTurb.ref()[i]
							* velocityProjectionNeighbour
							+ SNeighbour
									* (rhobState
											- surfaceNeighbourSide.rhobTurb.ref()[i]))
							* mesh.surfaces()[i].N();
		}
		else if (SNeighbour <= 0)
		{
			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceNeighbourSide.density[k].ref()[i];

			const vector momentumState = surfaceNeighbourSide.momentum.ref()[i];
			const scalar totalEnergyState =
					surfaceNeighbourSide.totalEnergy.ref()[i];
			const scalar rhokState = surfaceNeighbourSide.rhokTurb.ref()[i];
			const scalar rhoEpsState = surfaceNeighbourSide.rhoepsTurb.ref()[i];
			const vector rhoaState = surfaceNeighbourSide.rhoaTurb.ref()[i];
			const scalar rhobState = surfaceNeighbourSide.rhobTurb.ref()[i];

			velocityState = surfaceNeighbourSide.velocity.ref()[i];
			pressureState = surfaceNeighbourSide.pressure.ref()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].ref_r()[i] = velocityProjectionNeighbour
						* mesh.surfaces()[i].N() * densityState[k];

			numFluxes.momentum.ref_r()[i] = (momentumState
					* velocityProjectionNeighbour
					+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh.surfaces()[i].N())
							* (pressureState + turbulentPressureNeighbour))
					* mesh.surfaces()[i].N();

			numFluxes.totalEnergy.ref_r()[i] = velocityProjectionNeighbour
					* (totalEnergyState + pressureState
							+ turbulentPressureNeighbour)
					* mesh.surfaces()[i].N();

			numFluxes.rhokTurb.ref_r()[i] = velocityProjectionNeighbour
					* rhokState * mesh.surfaces()[i].N();
			numFluxes.rhoepsTurb.ref_r()[i] = velocityProjectionNeighbour
					* mesh.surfaces()[i].N() * rhoEpsState;
			numFluxes.rhoaTurb.ref_r()[i] = rhoaState
					* velocityProjectionNeighbour * mesh.surfaces()[i].N();
			numFluxes.rhobTurb.ref_r()[i] = velocityProjectionNeighbour
					* mesh.surfaces()[i].N() * rhobState;
		}
		else
		{
			errValue = errorsEnum::RiemannSolverError;
			throw exception("Fatal error in Riemann solver.", errValue);
		}

		starValues.c[0].ref_r()[i] = 0;
		for (std::size_t k = 1; k < densityState.size(); ++k)
		{
			starValues.c[k].ref_r()[i] = densityState[k]
					/ surfaceNeighbourSide.phaseThermodynamics->Mv()[k - 1];

			starValues.c[0].ref_r()[i] += starValues.c[k].ref()[i];
		}
		starValues.rho.ref_r()[i] = densityState[0];
		starValues.v.ref_r()[i] = velocityState;
		starValues.p.ref_r()[i] = pressureState;
		starValues.a.ref_r()[i] = aState;
		starValues.b.ref_r()[i] = bState;
	}

	return std::make_tuple(numFluxes, starValues);
}

