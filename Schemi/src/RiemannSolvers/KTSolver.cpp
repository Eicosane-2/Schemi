/*
 * KTSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "KTSolver.hpp"

#include "vector.hpp"

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
				densityState[k] = (surfaceNeighbourSide.density[k].cval()[i]
						* (SNeighbour - velocityProjectionNeighbour)
						- surfaceOwnerSide.density[k].cval()[i]
								* (SOwner - velocityProjectionOwner))
						* denominator;

			const vector momentumState =
					(surfaceNeighbourSide.momentum.cval()[i]
							* (SNeighbour - velocityProjectionNeighbour)
							- surfaceOwnerSide.momentum.cval()[i]
									* (SOwner - velocityProjectionOwner)
							- mesh_.surfaces()[i].N()
									* ((surfaceNeighbourSide.pressure.cval()[i]
											+ turbulentPressureNeighbour)
											- (surfaceOwnerSide.pressure.cval()[i]
													+ turbulentPressureOwner)))
							* denominator;
			const scalar totalEnergyState =
					(surfaceNeighbourSide.totalEnergy.cval()[i]
							* (SNeighbour - velocityProjectionNeighbour)
							- surfaceOwnerSide.totalEnergy.cval()[i]
									* (SOwner - velocityProjectionOwner)
							- ((surfaceNeighbourSide.pressure.cval()[i]
									+ turbulentPressureNeighbour)
									* velocityProjectionNeighbour
									- (surfaceOwnerSide.pressure.cval()[i]
											+ turbulentPressureOwner)
											* velocityProjectionOwner))
							* denominator;
			const scalar rhokState = (surfaceNeighbourSide.rhokTurb.cval()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhokTurb.cval()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;
			const vector rhoaState = (surfaceNeighbourSide.rhoaTurb.cval()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhoaTurb.cval()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;
			const scalar rhobState = (surfaceNeighbourSide.rhobTurb.cval()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhobTurb.cval()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;

			velocityState = momentumState / densityState[0];
			pressureState = pressureStar(
					*(surfaceOwnerSide.phaseThermodynamics), densityState,
					momentumState, totalEnergyState, rhokState);
			aState = rhoaState / densityState[0];
			bState = rhobState / densityState[0];

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].val()[i] =
						(SNeighbour * velocityProjectionOwner
								* surfaceOwnerSide.density[k].cval()[i]
								- SOwner * velocityProjectionNeighbour
										* surfaceNeighbourSide.density[k].cval()[i]
								+ SNeighbour * SOwner
										* (surfaceNeighbourSide.density[k].cval()[i]
												- surfaceOwnerSide.density[k].cval()[i]))
								* denominator * mesh_.surfaces()[i].N();

			numFluxes.momentum.val()[i] =
					(SNeighbour
							* (velocityProjectionOwner
									* surfaceOwnerSide.momentum.cval()[i]
									+ mesh_.surfaces()[i].N()
											* (surfaceOwnerSide.pressure.cval()[i]
													+ turbulentPressureOwner))
							- SOwner
									* (velocityProjectionNeighbour
											* surfaceNeighbourSide.momentum.cval()[i]
											+ mesh_.surfaces()[i].N()
													* (surfaceNeighbourSide.pressure.cval()[i]
															+ turbulentPressureNeighbour))
							+ SNeighbour * SOwner
									* (surfaceNeighbourSide.momentum.cval()[i]
											- surfaceOwnerSide.momentum.cval()[i]))
							* denominator * mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.val()[i] =
					(SNeighbour
							* (velocityProjectionOwner
									* (surfaceOwnerSide.totalEnergy.cval()[i]
											+ surfaceOwnerSide.pressure.cval()[i]
											+ turbulentPressureOwner))
							- SOwner
									* (velocityProjectionNeighbour
											* (surfaceNeighbourSide.totalEnergy.cval()[i]
													+ surfaceNeighbourSide.pressure.cval()[i]
													+ turbulentPressureNeighbour))
							+ SNeighbour * SOwner
									* (surfaceNeighbourSide.totalEnergy.cval()[i]
											- surfaceOwnerSide.totalEnergy.cval()[i]))
							* denominator * mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.val()[i] = (SNeighbour * velocityProjectionOwner
					* surfaceOwnerSide.rhokTurb.cval()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhokTurb.cval()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhokTurb.cval()[i]
									- surfaceOwnerSide.rhokTurb.cval()[i]))
					* denominator * mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.val()[i] = (SNeighbour
					* velocityProjectionOwner
					* surfaceOwnerSide.rhoepsTurb.cval()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhoepsTurb.cval()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhoepsTurb.cval()[i]
									- surfaceOwnerSide.rhoepsTurb.cval()[i]))
					* denominator * mesh_.surfaces()[i].N();
			numFluxes.rhoaTurb.val()[i] = (SNeighbour * velocityProjectionOwner
					* surfaceOwnerSide.rhoaTurb.cval()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhoaTurb.cval()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhoaTurb.cval()[i]
									- surfaceOwnerSide.rhoaTurb.cval()[i]))
					* denominator * mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.val()[i] = (SNeighbour * velocityProjectionOwner
					* surfaceOwnerSide.rhobTurb.cval()[i]
					- SOwner * velocityProjectionNeighbour
							* surfaceNeighbourSide.rhobTurb.cval()[i]
					+ SNeighbour * SOwner
							* (surfaceNeighbourSide.rhobTurb.cval()[i]
									- surfaceOwnerSide.rhobTurb.cval()[i]))
					* denominator * mesh_.surfaces()[i].N();
		}
		else
			[[unlikely]]
			throw exception("Fatal error in Riemann solver.", errValue =
					errors::RiemannSolverError);

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
