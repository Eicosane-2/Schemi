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
				densityState[k] = (surfaceNeighbourSide.density[k].ref()[i]
						* (SNeighbour - velocityProjectionNeighbour)
						- surfaceOwnerSide.density[k].ref()[i]
								* (SOwner - velocityProjectionOwner))
						* denominator;

			const vector momentumState = (surfaceNeighbourSide.momentum.ref()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.momentum.ref()[i]
							* (SOwner - velocityProjectionOwner)
					- mesh.surfaces()[i].N()
							* ((surfaceNeighbourSide.pressure.ref()[i]
									+ turbulentPressureNeighbour)
									- (surfaceOwnerSide.pressure.ref()[i]
											+ turbulentPressureOwner)))
					* denominator;
			const scalar totalEnergyState =
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
			const scalar rhokState = (surfaceNeighbourSide.rhokTurb.ref()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhokTurb.ref()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;
			const vector rhoaState = (surfaceNeighbourSide.rhoaTurb.ref()[i]
					* (SNeighbour - velocityProjectionNeighbour)
					- surfaceOwnerSide.rhoaTurb.ref()[i]
							* (SOwner - velocityProjectionOwner)) * denominator;
			const scalar rhobState = (surfaceNeighbourSide.rhobTurb.ref()[i]
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
		else
		{
			errValue = errors::RiemannSolverError;
			throw exception("Fatal error in Riemann solver.", errValue);
		}

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
