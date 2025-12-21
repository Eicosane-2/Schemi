/*
 * HLLCLMSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HLLCLMSolver.hpp"

#include "vector.hpp"
#include "tensor.hpp"

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

		scalar velocityProjectionRoe;
		scalar sqrSonicSpeedRoe;
		{
			const scalar sqrtRhoOwner { std::sqrt(
					surfaceOwnerSide.density[0].cval()[i]) };
			const scalar sqrtRhoNeighbour { std::sqrt(
					surfaceNeighbourSide.density[0].cval()[i]) };
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
					velocityProjectionOwner
							* (totalEnergyState + pressureState
									+ turbulentPressureOwner)
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

			for (std::size_t k = 0; k < densityStateO.size(); ++k)
				densityStateO[k] = surfaceOwnerSide.density[k].cval()[i]
						* factor * denominator;

			momentumStateO = (surfaceOwnerSide.momentum.cval()[i]
					- (surfaceOwnerSide.momentum.cval()[i]
							& mesh_.surfaces()[i].N()) * mesh_.surfaces()[i].N()
					+ surfaceOwnerSide.density[0].cval()[i] * Sc
							* mesh_.surfaces()[i].N()) * factor * denominator;

			totalEnergyStateO =
					factor * denominator
							* (surfaceOwnerSide.totalEnergy.cval()[i]
									+ (Sc - velocityProjectionOwner)
											* (surfaceOwnerSide.density[0].cval()[i]
													* Sc
													+ (surfaceOwnerSide.pressure.cval()[i]
															+ turbulentPressureOwner)
															/ factor));

			rhokStateO = surfaceOwnerSide.rhokTurb.cval()[i] * factor
					* denominator;
			rhoaStateO = surfaceOwnerSide.rhoaTurb.cval()[i] * factor
					* denominator;
			rhobStateO = surfaceOwnerSide.rhobTurb.cval()[i] * factor
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
				densityStateN[k] = surfaceNeighbourSide.density[k].cval()[i]
						* factor * denominator;

			momentumStateN = (surfaceNeighbourSide.momentum.cval()[i]
					- (surfaceNeighbourSide.momentum.cval()[i]
							& mesh_.surfaces()[i].N()) * mesh_.surfaces()[i].N()
					+ surfaceNeighbourSide.density[0].cval()[i] * Sc
							* mesh_.surfaces()[i].N()) * factor * denominator;

			totalEnergyStateN =
					factor * denominator
							* (surfaceNeighbourSide.totalEnergy.cval()[i]
									+ (Sc - velocityProjectionNeighbour)
											* (surfaceNeighbourSide.density[0].cval()[i]
													* Sc
													+ (surfaceNeighbourSide.pressure.cval()[i]
															+ turbulentPressureNeighbour)
															/ factor));

			rhokStateN = surfaceNeighbourSide.rhokTurb.cval()[i] * factor
					* denominator;
			rhoaStateN = surfaceNeighbourSide.rhoaTurb.cval()[i] * factor
					* denominator;
			rhobStateN = surfaceNeighbourSide.rhobTurb.cval()[i] * factor
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

		if (!(SOwner > 0) && !(SNeighbour < 0))
		{
			scalar rhoEpsStateO;
			scalar rhoEpsStateN;
			{
				const scalar factor { SOwner - velocityProjectionOwner };
				const scalar denominator { 1 / (SOwner - Sc) };

				for (std::size_t k = 0; k < densityStateO.size(); ++k)
					densityStateO[k] = surfaceOwnerSide.density[k].cval()[i]
							* factor * denominator;

				momentumStateO = (surfaceOwnerSide.momentum.cval()[i]
						- (surfaceOwnerSide.momentum.cval()[i]
								& mesh_.surfaces()[i].N())
								* mesh_.surfaces()[i].N()
						+ surfaceOwnerSide.density[0].cval()[i] * Sc
								* mesh_.surfaces()[i].N()) * factor
						* denominator;

				totalEnergyStateO =
						factor * denominator
								* (surfaceOwnerSide.totalEnergy.cval()[i]
										+ (Sc - velocityProjectionOwner)
												* (surfaceOwnerSide.density[0].cval()[i]
														* Sc
														+ (surfaceOwnerSide.pressure.cval()[i]
																+ turbulentPressureOwner)
																/ factor));

				rhokStateO = surfaceOwnerSide.rhokTurb.cval()[i] * factor
						* denominator;
				rhoEpsStateO = surfaceOwnerSide.rhoepsTurb.cval()[i] * factor
						* denominator;
				rhoaStateO = surfaceOwnerSide.rhoaTurb.cval()[i] * factor
						* denominator;
				rhobStateO = surfaceOwnerSide.rhobTurb.cval()[i] * factor
						* denominator;
			}
			{
				const scalar factor { SNeighbour - velocityProjectionNeighbour };
				const scalar denominator { 1 / (SNeighbour - Sc) };

				for (std::size_t k = 0; k < densityStateN.size(); ++k)
					densityStateN[k] = surfaceNeighbourSide.density[k].cval()[i]
							* factor * denominator;

				momentumStateN = (surfaceNeighbourSide.momentum.cval()[i]
						- (surfaceNeighbourSide.momentum.cval()[i]
								& mesh_.surfaces()[i].N())
								* mesh_.surfaces()[i].N()
						+ surfaceNeighbourSide.density[0].cval()[i] * Sc
								* mesh_.surfaces()[i].N()) * factor
						* denominator;

				totalEnergyStateN =
						factor * denominator
								* (surfaceNeighbourSide.totalEnergy.cval()[i]
										+ (Sc - velocityProjectionNeighbour)
												* (surfaceNeighbourSide.density[0].cval()[i]
														* Sc
														+ (surfaceNeighbourSide.pressure.cval()[i]
																+ turbulentPressureNeighbour)
																/ factor));

				rhokStateN = surfaceNeighbourSide.rhokTurb.cval()[i] * factor
						* denominator;
				rhoEpsStateN = surfaceNeighbourSide.rhoepsTurb.cval()[i]
						* factor * denominator;
				rhoaStateN = surfaceNeighbourSide.rhoaTurb.cval()[i] * factor
						* denominator;
				rhobStateN = surfaceNeighbourSide.rhobTurb.cval()[i] * factor
						* denominator;
			}

			const auto relMachNum = std::max(
					surfaceOwnerSide.velocity.cval()[i].mag()
							/ sonicSpeedOwner[i],
					surfaceNeighbourSide.velocity.cval()[i].mag()
							/ sonicSpeedNeighbour[i]) / minMachNumber;

			const auto phi = std::sin(
					std::min(static_cast<scalar>(1.), relMachNum) * Pi_number
							* 0.5);

			const scalar SOwner_corr { phi * SOwner };
			const scalar SNeighbour_corr { phi * SNeighbour };
			const scalar Sc_abs = std::abs(Sc);

			for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
				numFluxes.density[k].val()[i] =
						(0.5
								* (velocityProjectionOwner
										* surfaceOwnerSide.density[k].cval()[i]
										+ velocityProjectionNeighbour
												* surfaceNeighbourSide.density[k].cval()[i])
								+ 0.5
										* (SOwner_corr
												* (densityStateO[k]
														- surfaceOwnerSide.density[k].cval()[i])
												+ Sc_abs
														* (densityStateO[k]
																- densityStateN[k])
												+ SNeighbour_corr
														* (densityStateN[k]
																- surfaceNeighbourSide.density[k].cval()[i])))
								* mesh_.surfaces()[i].N();

			numFluxes.momentum.val()[i] =
					(0.5
							* (surfaceOwnerSide.momentum.cval()[i]
									* velocityProjectionOwner
									+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
											& mesh_.surfaces()[i].N())
											* (surfaceOwnerSide.pressure.cval()[i]
													+ turbulentPressureOwner)
									+ surfaceNeighbourSide.momentum.cval()[i]
											* velocityProjectionNeighbour
									+ (tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
											& mesh_.surfaces()[i].N())
											* (surfaceNeighbourSide.pressure.cval()[i]
													+ turbulentPressureNeighbour))
							+ 0.5
									* (SOwner_corr
											* (momentumStateO
													- surfaceOwnerSide.momentum.cval()[i])
											+ Sc_abs
													* (momentumStateO
															- momentumStateN)
											+ SNeighbour_corr
													* (momentumStateN
															- surfaceNeighbourSide.momentum.cval()[i])))
							* mesh_.surfaces()[i].N();

			numFluxes.totalEnergy.val()[i] =
					(0.5
							* (velocityProjectionOwner
									* (surfaceOwnerSide.totalEnergy.cval()[i]
											+ surfaceOwnerSide.pressure.cval()[i]
											+ turbulentPressureOwner)
									+ velocityProjectionNeighbour
											* (surfaceNeighbourSide.totalEnergy.cval()[i]
													+ surfaceNeighbourSide.pressure.cval()[i]
													+ turbulentPressureNeighbour))
							+ 0.5
									* (SOwner_corr
											* (totalEnergyStateO
													- surfaceOwnerSide.totalEnergy.cval()[i])
											+ Sc_abs
													* (totalEnergyStateO
															- totalEnergyStateN)
											+ SNeighbour_corr
													* (totalEnergyStateN
															- surfaceNeighbourSide.totalEnergy.cval()[i])))
							* mesh_.surfaces()[i].N();

			numFluxes.rhokTurb.val()[i] =
					(0.5
							* (velocityProjectionOwner
									* surfaceOwnerSide.rhokTurb.cval()[i]
									+ velocityProjectionNeighbour
											* surfaceNeighbourSide.rhokTurb.cval()[i])
							+ 0.5
									* (SOwner_corr
											* (rhokStateO
													- surfaceOwnerSide.rhokTurb.cval()[i])
											+ Sc_abs * (rhokStateO - rhokStateN)
											+ SNeighbour_corr
													* (rhokStateN
															- surfaceNeighbourSide.rhokTurb.cval()[i])))
							* mesh_.surfaces()[i].N();
			numFluxes.rhoepsTurb.val()[i] =
					(0.5
							* (velocityProjectionOwner
									* surfaceOwnerSide.rhoepsTurb.cval()[i]
									+ velocityProjectionNeighbour
											* surfaceNeighbourSide.rhoepsTurb.cval()[i])
							+ 0.5
									* (SOwner_corr
											* (rhoEpsStateO
													- surfaceOwnerSide.rhoepsTurb.cval()[i])
											+ Sc_abs
													* (rhoEpsStateO
															- rhoEpsStateN)
											+ SNeighbour_corr
													* (rhoEpsStateN
															- surfaceNeighbourSide.rhoepsTurb.cval()[i])))
							* mesh_.surfaces()[i].N();
			numFluxes.rhoaTurb.val()[i] =
					(0.5
							* (velocityProjectionOwner
									* surfaceOwnerSide.rhoaTurb.cval()[i]
									+ velocityProjectionNeighbour
											* surfaceNeighbourSide.rhoaTurb.cval()[i])
							+ 0.5
									* (SOwner_corr
											* (rhoaStateO
													- surfaceOwnerSide.rhoaTurb.cval()[i])
											+ Sc_abs * (rhoaStateO - rhoaStateN)
											+ SNeighbour_corr
													* (rhoaStateN
															- surfaceNeighbourSide.rhoaTurb.cval()[i])))
							* mesh_.surfaces()[i].N();
			numFluxes.rhobTurb.val()[i] =
					(0.5
							* (velocityProjectionOwner
									* surfaceOwnerSide.rhobTurb.cval()[i]
									+ velocityProjectionNeighbour
											* surfaceNeighbourSide.rhobTurb.cval()[i])
							+ 0.5
									* (SOwner_corr
											* (rhobStateO
													- surfaceOwnerSide.rhobTurb.cval()[i])
											+ Sc_abs * (rhobStateO - rhobStateN)
											+ SNeighbour_corr
													* (rhobStateN
															- surfaceNeighbourSide.rhobTurb.cval()[i])))
							* mesh_.surfaces()[i].N();
		}

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
