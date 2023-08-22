/*
 * HLLCLMSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HLLCLMSolver.hpp"

#include "dyadicProduct.hpp"
#include "tensorVectorDotProduct.hpp"
#include "vectorVectorDotProduct.hpp"

std::tuple<schemi::conservativeFlows, schemi::starFields> schemi::HLLCLMSolver::calculateFlows(
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

		scalar velocityProjectionRoe;
		scalar sqrSonicSpeedRoe;
		{
			const scalar sqrtRhoOwner { std::sqrt(
					surfaceOwnerSide.density[0].ref()[i]) };
			const scalar sqrtRhoNeighbour { std::sqrt(
					surfaceNeighbourSide.density[0].ref()[i]) };
			velocityProjectionRoe = (velocityProjectionOwner * sqrtRhoOwner
					+ velocityProjectionNeighbour * sqrtRhoNeighbour)
					/ (sqrtRhoOwner + sqrtRhoNeighbour);

			sqrSonicSpeedRoe = (sonicSpeedOwner[i] * sonicSpeedOwner[i]
					* sqrtRhoOwner
					+ sonicSpeedNeighbour[i] * sonicSpeedNeighbour[i]
							* sqrtRhoNeighbour)
					/ (sqrtRhoOwner + sqrtRhoNeighbour)
					+ 0.5 * (sqrtRhoOwner * sqrtRhoNeighbour)
							/ ((sqrtRhoOwner + sqrtRhoNeighbour)
									* (sqrtRhoOwner + sqrtRhoNeighbour))
							* (velocityProjectionNeighbour
									- velocityProjectionOwner)
							* (velocityProjectionNeighbour
									- velocityProjectionOwner);
		}

		const scalar SOwner { std::min(
				velocityProjectionOwner - sonicSpeedOwner[i],
				velocityProjectionRoe - std::sqrt(sqrSonicSpeedRoe)) };
		const scalar SNeighbour { std::max(
				velocityProjectionNeighbour + sonicSpeedNeighbour[i],
				velocityProjectionRoe + std::sqrt(sqrSonicSpeedRoe)) };

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

		std::valarray<scalar> densityState(numFluxes.density.size());
		vector velocityState;
		scalar pressureState;
		vector aState;
		scalar bState;

		std::valarray<scalar> densityStateO(numFluxes.density.size());
		vector momentumStateO;
		scalar totalEnergyStateO;
		scalar rhokStateO;
		vector rhoaStateO;
		scalar rhobStateO;

		std::valarray<scalar> densityStateN(numFluxes.density.size());
		vector momentumStateN;
		scalar totalEnergyStateN;
		scalar rhokStateN;
		vector rhoaStateN;
		scalar rhobStateN;

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
					velocityProjectionOwner
							* (totalEnergyState + pressureState
									+ turbulentPressureOwner)
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

			for (std::size_t k = 0; k < densityStateO.size(); ++k)
				densityStateO[k] = surfaceOwnerSide.density[k].ref()[i] * factor
						* denominator;

			momentumStateO = (surfaceOwnerSide.momentum.ref()[i]
					- (surfaceOwnerSide.momentum.ref()[i]
							& mesh.surfaces()[i].N()) * mesh.surfaces()[i].N()
					+ surfaceOwnerSide.density[0].ref()[i] * Sc
							* mesh.surfaces()[i].N()) * factor * denominator;

			totalEnergyStateO =
					factor * denominator
							* (surfaceOwnerSide.totalEnergy.ref()[i]
									+ (Sc - velocityProjectionOwner)
											* (surfaceOwnerSide.density[0].ref()[i]
													* Sc
													+ (surfaceOwnerSide.pressure.ref()[i]
															+ turbulentPressureOwner)
															/ factor));

			rhokStateO = surfaceOwnerSide.rhokTurb.ref()[i] * factor
					* denominator;
			rhoaStateO = surfaceOwnerSide.rhoaTurb.ref()[i] * factor
					* denominator;
			rhobStateO = surfaceOwnerSide.rhobTurb.ref()[i] * factor
					* denominator;

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = densityStateO[k];

			const vector momentumState = momentumStateO;
			const scalar totalEnergyState = totalEnergyStateO;
			const scalar rhokState = rhokStateO;
			const vector rhoaState = rhoaStateO;
			const scalar rhobState = rhobStateO;
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

			for (std::size_t k = 0; k < densityStateN.size(); ++k)
				densityStateN[k] = surfaceNeighbourSide.density[k].ref()[i]
						* factor * denominator;

			momentumStateN = (surfaceNeighbourSide.momentum.ref()[i]
					- (surfaceNeighbourSide.momentum.ref()[i]
							& mesh.surfaces()[i].N()) * mesh.surfaces()[i].N()
					+ surfaceNeighbourSide.density[0].ref()[i] * Sc
							* mesh.surfaces()[i].N()) * factor * denominator;

			totalEnergyStateN =
					factor * denominator
							* (surfaceNeighbourSide.totalEnergy.ref()[i]
									+ (Sc - velocityProjectionNeighbour)
											* (surfaceNeighbourSide.density[0].ref()[i]
													* Sc
													+ (surfaceNeighbourSide.pressure.ref()[i]
															+ turbulentPressureNeighbour)
															/ factor));

			rhokStateN = surfaceNeighbourSide.rhokTurb.ref()[i] * factor
					* denominator;
			rhoaStateN = surfaceNeighbourSide.rhoaTurb.ref()[i] * factor
					* denominator;
			rhobStateN = surfaceNeighbourSide.rhobTurb.ref()[i] * factor
					* denominator;

			for (std::size_t k = 0; k < densityState.size(); ++k)
				densityState[k] = densityStateN[k];

			const vector momentumState = momentumStateN;
			const scalar totalEnergyState = totalEnergyStateN;
			const scalar rhokState = rhokStateN;
			const vector rhoaState = rhoaStateN;
			const scalar rhobState = rhobStateN;

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
			errValue = errors::RiemannSolverError;
			throw exception("Fatal error in Riemann solver.", errValue);
		}

		if (!(SOwner > 0) && !(SNeighbour < 0))
		{
			scalar rhoEpsStateO;
			scalar rhoEpsStateN;
			{
				const scalar factor { SOwner - velocityProjectionOwner };
				const scalar denominator { 1 / (SOwner - Sc) };

				for (std::size_t k = 0; k < densityStateO.size(); ++k)
					densityStateO[k] = surfaceOwnerSide.density[k].ref()[i]
							* factor * denominator;

				momentumStateO = (surfaceOwnerSide.momentum.ref()[i]
						- (surfaceOwnerSide.momentum.ref()[i]
								& mesh.surfaces()[i].N())
								* mesh.surfaces()[i].N()
						+ surfaceOwnerSide.density[0].ref()[i] * Sc
								* mesh.surfaces()[i].N()) * factor
						* denominator;

				totalEnergyStateO =
						factor * denominator
								* (surfaceOwnerSide.totalEnergy.ref()[i]
										+ (Sc - velocityProjectionOwner)
												* (surfaceOwnerSide.density[0].ref()[i]
														* Sc
														+ (surfaceOwnerSide.pressure.ref()[i]
																+ turbulentPressureOwner)
																/ factor));

				rhokStateO = surfaceOwnerSide.rhokTurb.ref()[i] * factor
						* denominator;
				rhoEpsStateO = surfaceOwnerSide.rhoepsTurb.ref()[i] * factor
						* denominator;
				rhoaStateO = surfaceOwnerSide.rhoaTurb.ref()[i] * factor
						* denominator;
				rhobStateO = surfaceOwnerSide.rhobTurb.ref()[i] * factor
						* denominator;
			}
			{
				const scalar factor { SNeighbour - velocityProjectionNeighbour };
				const scalar denominator { 1 / (SNeighbour - Sc) };

				for (std::size_t k = 0; k < densityStateN.size(); ++k)
					densityStateN[k] = surfaceNeighbourSide.density[k].ref()[i]
							* factor * denominator;

				momentumStateN = (surfaceNeighbourSide.momentum.ref()[i]
						- (surfaceNeighbourSide.momentum.ref()[i]
								& mesh.surfaces()[i].N())
								* mesh.surfaces()[i].N()
						+ surfaceNeighbourSide.density[0].ref()[i] * Sc
								* mesh.surfaces()[i].N()) * factor
						* denominator;

				totalEnergyStateN =
						factor * denominator
								* (surfaceNeighbourSide.totalEnergy.ref()[i]
										+ (Sc - velocityProjectionNeighbour)
												* (surfaceNeighbourSide.density[0].ref()[i]
														* Sc
														+ (surfaceNeighbourSide.pressure.ref()[i]
																+ turbulentPressureNeighbour)
																/ factor));

				rhokStateN = surfaceNeighbourSide.rhokTurb.ref()[i] * factor
						* denominator;
				rhoEpsStateN = surfaceNeighbourSide.rhoepsTurb.ref()[i] * factor
						* denominator;
				rhoaStateN = surfaceNeighbourSide.rhoaTurb.ref()[i] * factor
						* denominator;
				rhobStateN = surfaceNeighbourSide.rhobTurb.ref()[i] * factor
						* denominator;
			}

			const auto relMachNum = std::max(
					surfaceOwnerSide.velocity.ref()[i].mag()
							/ sonicSpeedOwner[i],
					surfaceNeighbourSide.velocity.ref()[i].mag()
							/ sonicSpeedNeighbour[i]) / minMachNumber;

			const auto phi = std::sin(
					std::min(static_cast<scalar>(1.), relMachNum) * Pi_number
							* 0.5);

			const scalar SOwner_corr { phi * SOwner };
			const scalar SNeighbour_corr { phi * SNeighbour };
			const scalar Sc_abs = std::abs(Sc);

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].ref_r()[i] =
						(0.5
								* (velocityProjectionOwner
										* surfaceOwnerSide.density[k].ref()[i]
										+ velocityProjectionNeighbour
												* surfaceNeighbourSide.density[k].ref()[i])
								+ 0.5
										* (SOwner_corr
												* (densityStateO[k]
														- surfaceOwnerSide.density[k].ref()[i])
												+ Sc_abs
														* (densityStateO[k]
																- densityStateN[k])
												+ SNeighbour_corr
														* (densityStateN[k]
																- surfaceNeighbourSide.density[k].ref()[i])))
								* mesh.surfaces()[i].N();

			numFluxes.momentum.ref_r()[i] =
					(0.5
							* (surfaceOwnerSide.momentum.ref()[i]
									* velocityProjectionOwner
									+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
											& mesh.surfaces()[i].N())
											* (surfaceOwnerSide.pressure.ref()[i]
													+ turbulentPressureOwner)
									+ surfaceNeighbourSide.momentum.ref()[i]
											* velocityProjectionNeighbour
									+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
											& mesh.surfaces()[i].N())
											* (surfaceNeighbourSide.pressure.ref()[i]
													+ turbulentPressureNeighbour))
							+ 0.5
									* (SOwner_corr
											* (momentumStateO
													- surfaceOwnerSide.momentum.ref()[i])
											+ Sc_abs
													* (momentumStateO
															- momentumStateN)
											+ SNeighbour_corr
													* (momentumStateN
															- surfaceNeighbourSide.momentum.ref()[i])))
							* mesh.surfaces()[i].N();

			numFluxes.totalEnergy.ref_r()[i] =
					(0.5
							* (velocityProjectionOwner
									* (surfaceOwnerSide.totalEnergy.ref()[i]
											+ surfaceOwnerSide.pressure.ref()[i]
											+ turbulentPressureOwner)
									+ velocityProjectionNeighbour
											* (surfaceNeighbourSide.totalEnergy.ref()[i]
													+ surfaceNeighbourSide.pressure.ref()[i]
													+ turbulentPressureNeighbour))
							+ 0.5
									* (SOwner_corr
											* (totalEnergyStateO
													- surfaceOwnerSide.totalEnergy.ref()[i])
											+ Sc_abs
													* (totalEnergyStateO
															- totalEnergyStateN)
											+ SNeighbour_corr
													* (totalEnergyStateN
															- surfaceNeighbourSide.totalEnergy.ref()[i])))
							* mesh.surfaces()[i].N();

			numFluxes.rhokTurb.ref_r()[i] =
					(0.5
							* (velocityProjectionOwner
									* surfaceOwnerSide.rhokTurb.ref()[i]
									+ velocityProjectionNeighbour
											* surfaceNeighbourSide.rhokTurb.ref()[i])
							+ 0.5
									* (SOwner_corr
											* (rhokStateO
													- surfaceOwnerSide.rhokTurb.ref()[i])
											+ Sc_abs * (rhokStateO - rhokStateN)
											+ SNeighbour_corr
													* (rhokStateN
															- surfaceNeighbourSide.rhokTurb.ref()[i])))
							* mesh.surfaces()[i].N();
			numFluxes.rhoepsTurb.ref_r()[i] =
					(0.5
							* (velocityProjectionOwner
									* surfaceOwnerSide.rhoepsTurb.ref()[i]
									+ velocityProjectionNeighbour
											* surfaceNeighbourSide.rhoepsTurb.ref()[i])
							+ 0.5
									* (SOwner_corr
											* (rhoEpsStateO
													- surfaceOwnerSide.rhoepsTurb.ref()[i])
											+ Sc_abs
													* (rhoEpsStateO
															- rhoEpsStateN)
											+ SNeighbour_corr
													* (rhoEpsStateN
															- surfaceNeighbourSide.rhoepsTurb.ref()[i])))
							* mesh.surfaces()[i].N();
			numFluxes.rhoaTurb.ref_r()[i] =
					(0.5
							* (velocityProjectionOwner
									* surfaceOwnerSide.rhoaTurb.ref()[i]
									+ velocityProjectionNeighbour
											* surfaceNeighbourSide.rhoaTurb.ref()[i])
							+ 0.5
									* (SOwner_corr
											* (rhoaStateO
													- surfaceOwnerSide.rhoaTurb.ref()[i])
											+ Sc_abs * (rhoaStateO - rhoaStateN)
											+ SNeighbour_corr
													* (rhoaStateN
															- surfaceNeighbourSide.rhoaTurb.ref()[i])))
							* mesh.surfaces()[i].N();
			numFluxes.rhobTurb.ref_r()[i] =
					(0.5
							* (velocityProjectionOwner
									* surfaceOwnerSide.rhobTurb.ref()[i]
									+ velocityProjectionNeighbour
											* surfaceNeighbourSide.rhobTurb.ref()[i])
							+ 0.5
									* (SOwner_corr
											* (rhobStateO
													- surfaceOwnerSide.rhobTurb.ref()[i])
											+ Sc_abs * (rhobStateO - rhobStateN)
											+ SNeighbour_corr
													* (rhobStateN
															- surfaceNeighbourSide.rhobTurb.ref()[i])))
							* mesh.surfaces()[i].N();
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
