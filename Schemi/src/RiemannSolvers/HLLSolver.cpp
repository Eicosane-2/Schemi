/*
 * HLLSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HLLSolver.hpp"

#include "dyadicProduct.hpp"
#include "tensorVectorDotProduct.hpp"
#include "vectorVectorDotProduct.hpp"

std::tuple<schemi::conservativeFlows, schemi::starFields> schemi::HLLSolver::calculateFlows(
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

		std::valarray<scalar> densityState(numFluxes.density.size());
		vector velocityState;
		scalar pressureState;
		vector aState;
		scalar bState;

		if (SOwner >= 0)
		{
			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceOwnerSide.density[k].ref()[i];

			const auto & momentumState = surfaceOwnerSide.momentum.ref()[i];
			const auto & totalEnergyState =
					surfaceOwnerSide.totalEnergy.ref()[i];
			const auto & rhokState = surfaceOwnerSide.rhokTurb.ref()[i];
			const auto & rhoEpsState = surfaceOwnerSide.rhoepsTurb.ref()[i];
			const auto & rhoaState = surfaceOwnerSide.rhoaTurb.ref()[i];
			const auto & rhobState = surfaceOwnerSide.rhobTurb.ref()[i];

			velocityState = surfaceOwnerSide.velocity.ref()[i];
			pressureState = surfaceOwnerSide.pressure.ref()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].ref_r()[i] = velocityProjectionOwner
						* densityState[k] * mesh.surfaces()[i].N();

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
					* rhoEpsState * mesh.surfaces()[i].N();
			numFluxes.rhoaTurb.ref_r()[i] = rhoaState * velocityProjectionOwner
					* mesh.surfaces()[i].N();
			numFluxes.rhobTurb.ref_r()[i] = velocityProjectionOwner * rhobState
					* mesh.surfaces()[i].N();
		}
		else if (SNeighbour > 0)
		{
			const scalar denominator { 1 / (SNeighbour - SOwner) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = (surfaceNeighbourSide.density[k].ref()[i]
						* (SNeighbour - velocityProjectionNeighbour)
						- surfaceOwnerSide.density[k].ref()[i]
								* (SOwner - velocityProjectionOwner))
						* denominator;

			const auto momentumState = (surfaceNeighbourSide.momentum.ref()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.momentum.ref()[i]
							* (SOwner - velocityProjectionOwner)
					- mesh.surfaces()[i].N()
							* ((surfaceNeighbourSide.pressure.ref()[i]
									+ turbulentPressureNeighbour)
									- (surfaceOwnerSide.pressure.ref()[i]
											+ turbulentPressureOwner)))
					* denominator;
			const auto totalEnergyState =
					(surfaceNeighbourSide.totalEnergy.ref()[i]
							* (SNeighbour - velocityProjectionNeighbour)
							- surfaceOwnerSide.totalEnergy.ref()[i]
									* (SOwner - velocityProjectionOwner)
							- ((surfaceNeighbourSide.pressure.ref()[i]
									+ turbulentPressureNeighbour)
									* velocityProjectionNeighbour
									- (surfaceOwnerSide.pressure.ref()[i]
											+ turbulentPressureOwner)
											* velocityProjectionOwner))
							* denominator;
			const auto rhokState = (surfaceNeighbourSide.rhokTurb.ref()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhokTurb.ref()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;
			const auto rhoaState = (surfaceNeighbourSide.rhoaTurb.ref()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhoaTurb.ref()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;
			const auto rhobState = (surfaceNeighbourSide.rhobTurb.ref()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhobTurb.ref()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceOwnerSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].ref_r()[i] = (SNeighbour
						* velocityProjectionOwner
						* surfaceOwnerSide.density[k].ref()[i]
						- SOwner * velocityProjectionNeighbour
								* surfaceNeighbourSide.density[k].ref()[i]
						+ SNeighbour * SOwner
								* (surfaceNeighbourSide.density[k].ref()[i]
										- surfaceOwnerSide.density[k].ref()[i]))
						* denominator * mesh.surfaces()[i].N();

			numFluxes.momentum.ref_r()[i] =
					(SNeighbour
							* (velocityProjectionOwner
									* surfaceOwnerSide.momentum.ref()[i]
									+ mesh.surfaces()[i].N()
											* (surfaceOwnerSide.pressure.ref()[i]
													+ turbulentPressureOwner))
							- SOwner
									* (velocityProjectionNeighbour
											* surfaceNeighbourSide.momentum.ref()[i]
											+ mesh.surfaces()[i].N()
													* (surfaceNeighbourSide.pressure.ref()[i]
															+ turbulentPressureNeighbour))
							+ SNeighbour * SOwner
									* (surfaceNeighbourSide.momentum.ref()[i]
											- surfaceOwnerSide.momentum.ref()[i]))
							* denominator * mesh.surfaces()[i].N();

			numFluxes.totalEnergy.ref_r()[i] =
					(SNeighbour
							* (velocityProjectionOwner
									* (surfaceOwnerSide.totalEnergy.ref()[i]
											+ surfaceOwnerSide.pressure.ref()[i]
											+ turbulentPressureOwner))
							- SOwner
									* (velocityProjectionNeighbour
											* (surfaceNeighbourSide.totalEnergy.ref()[i]
													+ surfaceNeighbourSide.pressure.ref()[i]
													+ turbulentPressureNeighbour))
							+ SNeighbour * SOwner
									* (surfaceNeighbourSide.totalEnergy.ref()[i]
											- surfaceOwnerSide.totalEnergy.ref()[i]))
							* denominator * mesh.surfaces()[i].N();

			numFluxes.rhokTurb.ref_r()[i] = (SNeighbour
					* velocityProjectionOwner
					* surfaceOwnerSide.rhokTurb.ref()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhokTurb.ref()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhokTurb.ref()[i]
									- surfaceOwnerSide.rhokTurb.ref()[i]))
					* denominator * mesh.surfaces()[i].N();
			numFluxes.rhoepsTurb.ref_r()[i] = (SNeighbour
					* velocityProjectionOwner
					* surfaceOwnerSide.rhoepsTurb.ref()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhoepsTurb.ref()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhoepsTurb.ref()[i]
									- surfaceOwnerSide.rhoepsTurb.ref()[i]))
					* denominator * mesh.surfaces()[i].N();
			numFluxes.rhoaTurb.ref_r()[i] = (SNeighbour
					* velocityProjectionOwner
					* surfaceOwnerSide.rhoaTurb.ref()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhoaTurb.ref()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhoaTurb.ref()[i]
									- surfaceOwnerSide.rhoaTurb.ref()[i]))
					* denominator * mesh.surfaces()[i].N();
			numFluxes.rhobTurb.ref_r()[i] = (SNeighbour
					* velocityProjectionOwner
					* surfaceOwnerSide.rhobTurb.ref()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhobTurb.ref()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhobTurb.ref()[i]
									- surfaceOwnerSide.rhobTurb.ref()[i]))
					* denominator * mesh.surfaces()[i].N();
		}
		else if (SNeighbour <= 0)
		{
			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = surfaceNeighbourSide.density[k].ref()[i];

			const auto & momentumState = surfaceNeighbourSide.momentum.ref()[i];
			const auto & totalEnergyState =
					surfaceNeighbourSide.totalEnergy.ref()[i];
			const auto & rhokState = surfaceNeighbourSide.rhokTurb.ref()[i];
			const auto & rhoEpsState = surfaceNeighbourSide.rhoepsTurb.ref()[i];
			const auto & rhoaState = surfaceNeighbourSide.rhoaTurb.ref()[i];
			const auto & rhobState = surfaceNeighbourSide.rhobTurb.ref()[i];

			velocityState = surfaceNeighbourSide.velocity.ref()[i];
			pressureState = surfaceNeighbourSide.pressure.ref()[i];
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].ref_r()[i] = velocityProjectionNeighbour
						* densityState[k] * mesh.surfaces()[i].N();

			numFluxes.momentum.ref_r()[i] = (momentumState
					* velocityProjectionNeighbour
					+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
							& mesh.surfaces()[i].N())
							* (pressureState + turbulentPressureNeighbour))
					* mesh.surfaces()[i].N();

			numFluxes.totalEnergy.ref_r()[i] = (velocityProjectionNeighbour
					* (totalEnergyState + pressureState
							+ turbulentPressureNeighbour))
					* mesh.surfaces()[i].N();

			numFluxes.rhokTurb.ref_r()[i] = velocityProjectionNeighbour
					* rhokState * mesh.surfaces()[i].N();
			numFluxes.rhoepsTurb.ref_r()[i] = velocityProjectionNeighbour
					* rhoEpsState * mesh.surfaces()[i].N();
			numFluxes.rhoaTurb.ref_r()[i] = rhoaState
					* velocityProjectionNeighbour * mesh.surfaces()[i].N();
			numFluxes.rhobTurb.ref_r()[i] = velocityProjectionNeighbour
					* rhobState * mesh.surfaces()[i].N();
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
					/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

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
