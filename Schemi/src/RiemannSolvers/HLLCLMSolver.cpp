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

		scalar velocityProjectionRoe;
		scalar sqrSonicSpeedRoe;
		{
			const scalar sqrtRhoOwner { std::sqrt(
					surfaceOwnerSide.density[0]()[i]) };
			const scalar sqrtRhoNeighbour { std::sqrt(
					surfaceNeighbourSide.density[0]()[i]) };
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
					velocityProjectionOwner
							* (totalEnergyState + pressureState
									+ turbulentPressureOwner)
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

			for (std::size_t k = 0; k < densityStateO.size(); ++k)
				densityStateO[k] = surfaceOwnerSide.density[k]()[i] * factor
						* denominator;

			momentumStateO = (surfaceOwnerSide.momentum()[i]
					- (surfaceOwnerSide.momentum()[i] & mesh_.surfaces()[i].N())
							* mesh_.surfaces()[i].N()
					+ surfaceOwnerSide.density[0]()[i] * Sc
							* mesh_.surfaces()[i].N()) * factor * denominator;

			totalEnergyStateO = factor * denominator
					* (surfaceOwnerSide.totalEnergy()[i]
							+ (Sc - velocityProjectionOwner)
									* (surfaceOwnerSide.density[0]()[i] * Sc
											+ (surfaceOwnerSide.pressure()[i]
													+ turbulentPressureOwner)
													/ factor));

			rhokStateO = surfaceOwnerSide.rhokTurb()[i] * factor * denominator;
			rhoaStateO = surfaceOwnerSide.rhoaTurb()[i] * factor * denominator;
			rhobStateO = surfaceOwnerSide.rhobTurb()[i] * factor * denominator;

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
				densityStateN[k] = surfaceNeighbourSide.density[k]()[i] * factor
						* denominator;

			momentumStateN = (surfaceNeighbourSide.momentum()[i]
					- (surfaceNeighbourSide.momentum()[i]
							& mesh_.surfaces()[i].N()) * mesh_.surfaces()[i].N()
					+ surfaceNeighbourSide.density[0]()[i] * Sc
							* mesh_.surfaces()[i].N()) * factor * denominator;

			totalEnergyStateN =
					factor * denominator
							* (surfaceNeighbourSide.totalEnergy()[i]
									+ (Sc - velocityProjectionNeighbour)
											* (surfaceNeighbourSide.density[0]()[i]
													* Sc
													+ (surfaceNeighbourSide.pressure()[i]
															+ turbulentPressureNeighbour)
															/ factor));

			rhokStateN = surfaceNeighbourSide.rhokTurb()[i] * factor
					* denominator;
			rhoaStateN = surfaceNeighbourSide.rhoaTurb()[i] * factor
					* denominator;
			rhobStateN = surfaceNeighbourSide.rhobTurb()[i] * factor
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
			[[unlikely]]
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
						densityStateO[k] = surfaceOwnerSide.density[k]()[i]
								* factor * denominator;

					momentumStateO = (surfaceOwnerSide.momentum()[i]
							- (surfaceOwnerSide.momentum()[i]
									& mesh_.surfaces()[i].N())
									* mesh_.surfaces()[i].N()
							+ surfaceOwnerSide.density[0]()[i] * Sc
									* mesh_.surfaces()[i].N()) * factor
							* denominator;

					totalEnergyStateO =
							factor * denominator
									* (surfaceOwnerSide.totalEnergy()[i]
											+ (Sc - velocityProjectionOwner)
													* (surfaceOwnerSide.density[0]()[i]
															* Sc
															+ (surfaceOwnerSide.pressure()[i]
																	+ turbulentPressureOwner)
																	/ factor));

					rhokStateO = surfaceOwnerSide.rhokTurb()[i] * factor
							* denominator;
					rhoEpsStateO = surfaceOwnerSide.rhoepsTurb()[i] * factor
							* denominator;
					rhoaStateO = surfaceOwnerSide.rhoaTurb()[i] * factor
							* denominator;
					rhobStateO = surfaceOwnerSide.rhobTurb()[i] * factor
							* denominator;
				}
				{
					const scalar factor { SNeighbour
							- velocityProjectionNeighbour };
					const scalar denominator { 1 / (SNeighbour - Sc) };

					for (std::size_t k = 0; k < densityStateN.size(); ++k)
						densityStateN[k] = surfaceNeighbourSide.density[k]()[i]
								* factor * denominator;

					momentumStateN = (surfaceNeighbourSide.momentum()[i]
							- (surfaceNeighbourSide.momentum()[i]
									& mesh_.surfaces()[i].N())
									* mesh_.surfaces()[i].N()
							+ surfaceNeighbourSide.density[0]()[i] * Sc
									* mesh_.surfaces()[i].N()) * factor
							* denominator;

					totalEnergyStateN =
							factor * denominator
									* (surfaceNeighbourSide.totalEnergy()[i]
											+ (Sc - velocityProjectionNeighbour)
													* (surfaceNeighbourSide.density[0]()[i]
															* Sc
															+ (surfaceNeighbourSide.pressure()[i]
																	+ turbulentPressureNeighbour)
																	/ factor));

					rhokStateN = surfaceNeighbourSide.rhokTurb()[i] * factor
							* denominator;
					rhoEpsStateN = surfaceNeighbourSide.rhoepsTurb()[i] * factor
							* denominator;
					rhoaStateN = surfaceNeighbourSide.rhoaTurb()[i] * factor
							* denominator;
					rhobStateN = surfaceNeighbourSide.rhobTurb()[i] * factor
							* denominator;
				}

				const auto relMachNum = std::max(
						surfaceOwnerSide.velocity()[i].mag()
								/ sonicSpeedOwner[i],
						surfaceNeighbourSide.velocity()[i].mag()
								/ sonicSpeedNeighbour[i]) / minMachNumber;

				const auto phi = std::sin(
						std::min(static_cast<scalar>(1.), relMachNum)
								* Pi_number * 0.5);

				const scalar SOwner_corr { phi * SOwner };
				const scalar SNeighbour_corr { phi * SNeighbour };
				const scalar Sc_abs = std::abs(Sc);

				for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
					numFluxes.density[k].r()[i] =
							(0.5
									* (velocityProjectionOwner
											* surfaceOwnerSide.density[k]()[i]
											+ velocityProjectionNeighbour
													* surfaceNeighbourSide.density[k]()[i])
									+ 0.5
											* (SOwner_corr
													* (densityStateO[k]
															- surfaceOwnerSide.density[k]()[i])
													+ Sc_abs
															* (densityStateO[k]
																	- densityStateN[k])
													+ SNeighbour_corr
															* (densityStateN[k]
																	- surfaceNeighbourSide.density[k]()[i])))
									* mesh_.surfaces()[i].N();

				numFluxes.momentum.r()[i] =
						(0.5
								* (surfaceOwnerSide.momentum()[i]
										* velocityProjectionOwner
										+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
												& mesh_.surfaces()[i].N())
												* (surfaceOwnerSide.pressure()[i]
														+ turbulentPressureOwner)
										+ surfaceNeighbourSide.momentum()[i]
												* velocityProjectionNeighbour
										+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
												& mesh_.surfaces()[i].N())
												* (surfaceNeighbourSide.pressure()[i]
														+ turbulentPressureNeighbour))
								+ 0.5
										* (SOwner_corr
												* (momentumStateO
														- surfaceOwnerSide.momentum()[i])
												+ Sc_abs
														* (momentumStateO
																- momentumStateN)
												+ SNeighbour_corr
														* (momentumStateN
																- surfaceNeighbourSide.momentum()[i])))
								* mesh_.surfaces()[i].N();

				numFluxes.totalEnergy.r()[i] =
						(0.5
								* (velocityProjectionOwner
										* (surfaceOwnerSide.totalEnergy()[i]
												+ surfaceOwnerSide.pressure()[i]
												+ turbulentPressureOwner)
										+ velocityProjectionNeighbour
												* (surfaceNeighbourSide.totalEnergy()[i]
														+ surfaceNeighbourSide.pressure()[i]
														+ turbulentPressureNeighbour))
								+ 0.5
										* (SOwner_corr
												* (totalEnergyStateO
														- surfaceOwnerSide.totalEnergy()[i])
												+ Sc_abs
														* (totalEnergyStateO
																- totalEnergyStateN)
												+ SNeighbour_corr
														* (totalEnergyStateN
																- surfaceNeighbourSide.totalEnergy()[i])))
								* mesh_.surfaces()[i].N();

				numFluxes.rhokTurb.r()[i] =
						(0.5
								* (velocityProjectionOwner
										* surfaceOwnerSide.rhokTurb()[i]
										+ velocityProjectionNeighbour
												* surfaceNeighbourSide.rhokTurb()[i])
								+ 0.5
										* (SOwner_corr
												* (rhokStateO
														- surfaceOwnerSide.rhokTurb()[i])
												+ Sc_abs
														* (rhokStateO
																- rhokStateN)
												+ SNeighbour_corr
														* (rhokStateN
																- surfaceNeighbourSide.rhokTurb()[i])))
								* mesh_.surfaces()[i].N();
				numFluxes.rhoepsTurb.r()[i] =
						(0.5
								* (velocityProjectionOwner
										* surfaceOwnerSide.rhoepsTurb()[i]
										+ velocityProjectionNeighbour
												* surfaceNeighbourSide.rhoepsTurb()[i])
								+ 0.5
										* (SOwner_corr
												* (rhoEpsStateO
														- surfaceOwnerSide.rhoepsTurb()[i])
												+ Sc_abs
														* (rhoEpsStateO
																- rhoEpsStateN)
												+ SNeighbour_corr
														* (rhoEpsStateN
																- surfaceNeighbourSide.rhoepsTurb()[i])))
								* mesh_.surfaces()[i].N();
				numFluxes.rhoaTurb.r()[i] =
						(0.5
								* (velocityProjectionOwner
										* surfaceOwnerSide.rhoaTurb()[i]
										+ velocityProjectionNeighbour
												* surfaceNeighbourSide.rhoaTurb()[i])
								+ 0.5
										* (SOwner_corr
												* (rhoaStateO
														- surfaceOwnerSide.rhoaTurb()[i])
												+ Sc_abs
														* (rhoaStateO
																- rhoaStateN)
												+ SNeighbour_corr
														* (rhoaStateN
																- surfaceNeighbourSide.rhoaTurb()[i])))
								* mesh_.surfaces()[i].N();
				numFluxes.rhobTurb.r()[i] =
						(0.5
								* (velocityProjectionOwner
										* surfaceOwnerSide.rhobTurb()[i]
										+ velocityProjectionNeighbour
												* surfaceNeighbourSide.rhobTurb()[i])
								+ 0.5
										* (SOwner_corr
												* (rhobStateO
														- surfaceOwnerSide.rhobTurb()[i])
												+ Sc_abs
														* (rhobStateO
																- rhobStateN)
												+ SNeighbour_corr
														* (rhobStateN
																- surfaceNeighbourSide.rhobTurb()[i])))
								* mesh_.surfaces()[i].N();
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
