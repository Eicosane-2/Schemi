/*
 * KTSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "KTSolver.hpp"

#include "dyadicProduct.hpp"
#include "vectorVectorDotProduct.hpp"

std::tuple<schemi::conservativeFlows, schemi::starFields> schemi::KTSolver::calculateFlows(
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

		scalar SOwner { std::min(velocityProjectionOwner - sonicSpeedOwner[i],
				velocityProjectionNeighbour - sonicSpeedNeighbour[i]) };
		scalar SNeighbour { std::max(
				velocityProjectionOwner + sonicSpeedOwner[i],
				velocityProjectionNeighbour + sonicSpeedNeighbour[i]) };

		const scalar SOwnerSNeighbourMag { std::max(std::abs(SOwner),
				std::abs(SNeighbour)) };

		SOwner = -SOwnerSNeighbourMag;
		SNeighbour = SOwnerSNeighbourMag;

		std::valarray<scalar> densityState(numFluxes.density.size());
		vector velocityState;
		scalar pressureState;
		vector aState;
		scalar bState;

		if ((SNeighbour >= 0) && (SOwner <= 0))
		{
			const scalar denominator { 1 / (SNeighbour - SOwner) };

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = (surfaceNeighbourSide.density[k]()[i]
						* (SNeighbour - velocityProjectionNeighbour)
						- surfaceOwnerSide.density[k]()[i]
								* (SOwner - velocityProjectionOwner))
						* denominator;

			const vector momentumState = (surfaceNeighbourSide.momentum()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.momentum()[i]
							* (SOwner - velocityProjectionOwner)
					- mesh_.surfaces()[i].N()
							* ((surfaceNeighbourSide.pressure()[i]
									+ turbulentPressureNeighbour)
									- (surfaceOwnerSide.pressure()[i]
											+ turbulentPressureOwner)))
					* denominator;
			const scalar totalEnergyState =
					(surfaceNeighbourSide.totalEnergy()[i]
							* (SNeighbour - velocityProjectionNeighbour)
							- surfaceOwnerSide.totalEnergy()[i]
									* (SOwner - velocityProjectionOwner)
							- ((surfaceNeighbourSide.pressure()[i]
									+ turbulentPressureNeighbour)
									* velocityProjectionNeighbour
									- (surfaceOwnerSide.pressure()[i]
											+ turbulentPressureOwner)
											* velocityProjectionOwner))
							* denominator;
			const scalar rhokState = (surfaceNeighbourSide.rhokTurb()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhokTurb()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;
			const vector rhoaState = (surfaceNeighbourSide.rhoaTurb()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhoaTurb()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;
			const scalar rhobState = (surfaceNeighbourSide.rhobTurb()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhobTurb()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceOwnerSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].r()[i] = (SNeighbour
						* velocityProjectionOwner
						* surfaceOwnerSide.density[k]()[i]
						- SOwner * velocityProjectionNeighbour
								* surfaceNeighbourSide.density[k]()[i]
						+ SNeighbour * SOwner
								* (surfaceNeighbourSide.density[k]()[i]
										- surfaceOwnerSide.density[k]()[i]))
						* denominator * mesh_.surfaces()[i].N();

			numFluxes.momentum.r()[i] =
					(SNeighbour
							* (velocityProjectionOwner
									* surfaceOwnerSide.momentum()[i]
									+ mesh_.surfaces()[i].N()
											* (surfaceOwnerSide.pressure()[i]
													+ turbulentPressureOwner))
							- SOwner
									* (velocityProjectionNeighbour
											* surfaceNeighbourSide.momentum()[i]
											+ mesh_.surfaces()[i].N()
													* (surfaceNeighbourSide.pressure()[i]
															+ turbulentPressureNeighbour))
							+ SNeighbour * SOwner
									* (surfaceNeighbourSide.momentum()[i]
											- surfaceOwnerSide.momentum()[i]))
							* denominator * mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.r()[i] = (SNeighbour
					* (velocityProjectionOwner
							* (surfaceOwnerSide.totalEnergy()[i]
									+ surfaceOwnerSide.pressure()[i]
									+ turbulentPressureOwner))
					- SOwner
							* (velocityProjectionNeighbour
									* (surfaceNeighbourSide.totalEnergy()[i]
											+ surfaceNeighbourSide.pressure()[i]
											+ turbulentPressureNeighbour))
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.totalEnergy()[i]
									- surfaceOwnerSide.totalEnergy()[i]))
					* denominator * mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.r()[i] = (SNeighbour * velocityProjectionOwner
					* surfaceOwnerSide.rhokTurb()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhokTurb()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhokTurb()[i]
									- surfaceOwnerSide.rhokTurb()[i]))
					* denominator * mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.r()[i] = (SNeighbour * velocityProjectionOwner
					* surfaceOwnerSide.rhoepsTurb()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhoepsTurb()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhoepsTurb()[i]
									- surfaceOwnerSide.rhoepsTurb()[i]))
					* denominator * mesh_.surfaces()[i].N();
			numFluxes.rhoaTurb.r()[i] = (SNeighbour * velocityProjectionOwner
					* surfaceOwnerSide.rhoaTurb()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhoaTurb()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhoaTurb()[i]
									- surfaceOwnerSide.rhoaTurb()[i]))
					* denominator * mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.r()[i] = (SNeighbour * velocityProjectionOwner
					* surfaceOwnerSide.rhobTurb()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhobTurb()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhobTurb()[i]
									- surfaceOwnerSide.rhobTurb()[i]))
					* denominator * mesh_.surfaces()[i].N();
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
