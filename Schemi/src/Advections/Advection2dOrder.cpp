/*
 * Advection2dOrder.cpp
 *
 *  Created on: 2024/11/20
 *      Author: Maxim Boldyrev
 */

#include "Advection2dOrder.hpp"

#include <chrono>

#include "divergence.hpp"
#include "gradient.hpp"
#include "TVDLimiter.hpp"

schemi::starFields schemi::Advection2dOrder(
		homogeneousPhase<cubicCell> & gasPhase, const abstractLimiter & limiter,
		const abstractFlowSolver & fsolver, std::pair<bool, vector> gravitation,
		const boundaryConditionValue & boundaryConditionValueCalc,
		scalar & timeForTVD, scalar & timeForHancock,
		scalar & timeForFlowCalculation, scalar & timeForTimeIntegration,
		const MPIHandler & parallelism)
{
	auto & mesh_ { gasPhase.pressure.meshRef() };

	/*Creating surface fields.*/
	homogeneousPhase<quadraticSurface> surfaceOwnerSide { bunchOfFields<
			quadraticSurface>(gasPhase),
			transportCoefficients<quadraticSurface>(mesh_),
			gasPhase.phaseThermodynamics, gasPhase.turbulence,
			gasPhase.transportModel };

	homogeneousPhase<quadraticSurface> surfaceNeighbourSide { bunchOfFields<
			quadraticSurface>(gasPhase),
			transportCoefficients<quadraticSurface>(mesh_),
			gasPhase.phaseThermodynamics, gasPhase.turbulence,
			gasPhase.transportModel };

	/*Creating limited gradients.*/
	std::vector<volumeField<vector>> concentrationTVDGradient {
			gasPhase.phaseThermodynamics->Mv().size(), volumeField<vector>(
					mesh_, vector(0)) };
	volumeField<tensor> velocityTVDGradient { mesh_, tensor(0) };
	volumeField<vector> pressureTVDGradient { mesh_, vector(0) };
	volumeField<vector> kTVDGradient { mesh_, vector(0) };
	volumeField<vector> epsilonTVDGradient { mesh_, vector(0) };
	volumeField<tensor> aTVDGradient { mesh_, tensor(0) };
	volumeField<vector> bTVDGradient { mesh_, vector(0) };

	/*TVD Reconstruction.*/
	const auto TVDStartTime = std::chrono::high_resolution_clock::now();
	{
		/*TVD limiters.*/
		for (std::size_t k = 0; k < concentrationTVDGradient.size(); ++k)
		{
			const auto concentrationGradient_k = grad(
					gasPhase.concentration.v[k + 1], boundaryConditionValueCalc,
					k + 1);

			concentrationTVDGradient[k] = TVDLimiter(concentrationGradient_k,
					gasPhase.concentration.v[k + 1], limiter,
					boundaryConditionValueCalc, parallelism, k + 1);
		}

		const auto velocityGradient = grad(gasPhase.velocity,
				boundaryConditionValueCalc);

		velocityTVDGradient = TVDLimiter(velocityGradient, gasPhase.velocity,
				limiter, boundaryConditionValueCalc, parallelism);

		const auto pressureGradient = grad(gasPhase.pressure,
				boundaryConditionValueCalc);

		pressureTVDGradient = TVDLimiter(pressureGradient, gasPhase.pressure,
				limiter, boundaryConditionValueCalc, parallelism);

		if (gasPhase.turbulence->turbulence())
		{
			const auto kGradient = grad(gasPhase.kTurb,
					boundaryConditionValueCalc);

			kTVDGradient = TVDLimiter(kGradient, gasPhase.kTurb, limiter,
					boundaryConditionValueCalc, parallelism);

			const auto epsilonGradient = grad(gasPhase.epsTurb,
					boundaryConditionValueCalc);

			epsilonTVDGradient = TVDLimiter(epsilonGradient, gasPhase.epsTurb,
					limiter, boundaryConditionValueCalc, parallelism);

			if (gasPhase.turbulence->aField())
			{
				const auto aGradient = grad(gasPhase.aTurb,
						boundaryConditionValueCalc);

				aTVDGradient = TVDLimiter(aGradient, gasPhase.aTurb, limiter,
						boundaryConditionValueCalc, parallelism);

				if (gasPhase.turbulence->bField())
				{
					const auto bGradient = grad(gasPhase.bTurb,
							boundaryConditionValueCalc);

					bTVDGradient = TVDLimiter(bGradient, gasPhase.bTurb,
							limiter, boundaryConditionValueCalc, parallelism);
				}
			}
		}

		for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		{
			const std::vector<std::size_t> & surfacesOfCell_i {
					mesh_.surfacesOfCells()[i] };

			for (std::size_t j = 0; j < surfacesOfCell_i.size(); ++j)
			{
				const std::size_t surfaceIndex { surfacesOfCell_i[j] };

				const vector deltaVector = mesh_.surfaces()[surfaceIndex].rC()
						- mesh_.cells()[i].rC();

				std::valarray<scalar> reconstructedConcentrationsValues(
						gasPhase.phaseThermodynamics->Mv().size());

				for (std::size_t k = 0;
						k < reconstructedConcentrationsValues.size(); ++k)
					reconstructedConcentrationsValues[k] =
							(concentrationTVDGradient[k].cval()[i] & deltaVector)
									+ gasPhase.concentration.v[k + 1].cval()[i];

				const vector reconstructedVelocityValue =
						(velocityTVDGradient.cval()[i] & deltaVector)
								+ gasPhase.velocity.cval()[i];

				const scalar reconstructedPressureValue {
						(pressureTVDGradient.cval()[i] & deltaVector)
								+ gasPhase.pressure.cval()[i] };

				scalar reconstructedkValue { 0 };
				scalar reconstructedEpsilon { 0 };
				vector reconstructedaValue { 0 };
				scalar reconstructedbValue { 0 };
				if (gasPhase.turbulence->turbulence())
				{
					reconstructedkValue = (kTVDGradient.cval()[i] & deltaVector)
							+ gasPhase.kTurb.cval()[i];

					reconstructedEpsilon = (epsilonTVDGradient.cval()[i]
							& deltaVector) + gasPhase.epsTurb.cval()[i];

					if (gasPhase.turbulence->aField())
					{
						reconstructedaValue = (aTVDGradient.cval()[i]
								& deltaVector) + gasPhase.aTurb.cval()[i];

						if (gasPhase.turbulence->bField())
							reconstructedbValue = (bTVDGradient.cval()[i]
									& deltaVector) + gasPhase.bTurb.cval()[i];
					}
				}

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					for (std::size_t k = 0;
							k < reconstructedConcentrationsValues.size(); ++k)
						surfaceOwnerSide.concentration.v[k + 1].val()[surfaceIndex] =
								reconstructedConcentrationsValues[k];

					surfaceOwnerSide.velocity.val()[surfaceIndex] =
							reconstructedVelocityValue;

					surfaceOwnerSide.pressure.val()[surfaceIndex] =
							reconstructedPressureValue;

					if (gasPhase.turbulence->turbulence())
					{
						surfaceOwnerSide.kTurb.val()[surfaceIndex] =
								reconstructedkValue;

						surfaceOwnerSide.epsTurb.val()[surfaceIndex] =
								reconstructedEpsilon;

						if (gasPhase.turbulence->aField())
						{
							surfaceOwnerSide.aTurb.val()[surfaceIndex] =
									reconstructedaValue;

							if (gasPhase.turbulence->bField())
								surfaceOwnerSide.bTurb.val()[surfaceIndex] =
										reconstructedbValue;
						}
					}
				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					for (std::size_t k = 0;
							k < reconstructedConcentrationsValues.size(); ++k)
						surfaceNeighbourSide.concentration.v[k + 1].val()[surfaceIndex] =
								reconstructedConcentrationsValues[k];

					surfaceNeighbourSide.velocity.val()[surfaceIndex] =
							reconstructedVelocityValue;

					surfaceNeighbourSide.pressure.val()[surfaceIndex] =
							reconstructedPressureValue;

					if (gasPhase.turbulence->turbulence())
					{
						surfaceNeighbourSide.kTurb.val()[surfaceIndex] =
								reconstructedkValue;

						surfaceNeighbourSide.epsTurb.val()[surfaceIndex] =
								reconstructedEpsilon;

						if (gasPhase.turbulence->aField())
						{
							surfaceNeighbourSide.aTurb.val()[surfaceIndex] =
									reconstructedaValue;

							if (gasPhase.turbulence->bField())
								surfaceNeighbourSide.bTurb.val()[surfaceIndex] =
										reconstructedbValue;
						}
					}
				}
				else
					[[unlikely]]
					throw exception("Couldn't choose side to add",
							errors::systemError);
			}
		}

		/*Recalculation of other quantities for owner side*/
		surfaceOwnerSide.concentration.v[0].val() = 0;
		surfaceOwnerSide.density[0].val() = 0;
		for (std::size_t k = 1; k < surfaceOwnerSide.concentration.v.size();
				++k)
		{
			surfaceOwnerSide.concentration.v[0] +=
					surfaceOwnerSide.concentration.v[k];

			surfaceOwnerSide.density[k].val() =
					(surfaceOwnerSide.concentration.v[k]
							* surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1]).cval();

			surfaceOwnerSide.density[0] += surfaceOwnerSide.density[k];
		}

		if (gasPhase.turbulence->turbulence())
		{
			surfaceOwnerSide.rhokTurb.val() = (surfaceOwnerSide.kTurb
					* surfaceOwnerSide.density[0]).cval();

			surfaceOwnerSide.rhoepsTurb.val() = (surfaceOwnerSide.epsTurb
					* surfaceOwnerSide.density[0]).cval();

			if (gasPhase.turbulence->aField())
			{
				surfaceOwnerSide.rhoaTurb.val() = (surfaceOwnerSide.aTurb
						* surfaceOwnerSide.density[0]).cval();

				if (gasPhase.turbulence->bField())
					surfaceOwnerSide.rhobTurb.val() = (surfaceOwnerSide.bTurb
							* surfaceOwnerSide.density[0]).cval();
			}
		}

		surfaceOwnerSide.internalEnergy.val() =
				surfaceOwnerSide.phaseThermodynamics->UvFromp(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.pressure.cval());

		surfaceOwnerSide.temperature.val() =
				surfaceOwnerSide.phaseThermodynamics->TFromUv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.internalEnergy.cval());

		surfaceOwnerSide.momentum.val() = (surfaceOwnerSide.velocity
				* surfaceOwnerSide.density[0]).cval();

		{
			const auto v2 = surfaceOwnerSide.velocity
					& surfaceOwnerSide.velocity;

			surfaceOwnerSide.totalEnergy.val() =
					(surfaceOwnerSide.internalEnergy
							+ surfaceOwnerSide.density[0] * v2 * 0.5
							+ surfaceOwnerSide.rhokTurb).cval();
		}

		surfaceOwnerSide.HelmholtzEnergy.val() =
				surfaceOwnerSide.phaseThermodynamics->Fv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature.cval());

		surfaceOwnerSide.entropy.val() =
				surfaceOwnerSide.phaseThermodynamics->Sv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature.cval());

		/*Recalculation of other quantities for neighbour side*/
		const std::size_t nonExistentCell = mesh_.nonexistCell();
		for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
			if (mesh_.surfaceNeighbour()[i] != nonExistentCell)
			{
				surfaceNeighbourSide.concentration.v[0].val()[i] = 0;
				surfaceNeighbourSide.density[0].val()[i] = 0;
				for (std::size_t k = 1;
						k < surfaceNeighbourSide.concentration.v.size(); ++k)
				{
					surfaceNeighbourSide.concentration.v[0].val()[i] +=
							surfaceNeighbourSide.concentration.v[k].cval()[i];

					surfaceNeighbourSide.density[k].val()[i] =
							surfaceNeighbourSide.concentration.v[k].cval()[i]
									* gasPhase.phaseThermodynamics->Mv()[k - 1];

					surfaceNeighbourSide.density[0].val()[i] +=
							surfaceNeighbourSide.density[k].cval()[i];
				}

				surfaceNeighbourSide.momentum.val()[i] =
						surfaceNeighbourSide.velocity.cval()[i]
								* surfaceNeighbourSide.density[0].cval()[i];

				if (gasPhase.turbulence->turbulence())
				{
					surfaceNeighbourSide.rhokTurb.val()[i] =
							surfaceNeighbourSide.kTurb.cval()[i]
									* surfaceNeighbourSide.density[0].cval()[i];

					surfaceNeighbourSide.rhoepsTurb.val()[i] =
							surfaceNeighbourSide.epsTurb.cval()[i]
									* surfaceNeighbourSide.density[0].cval()[i];

					if (gasPhase.turbulence->aField())
					{
						surfaceNeighbourSide.rhoaTurb.val()[i] =
								surfaceNeighbourSide.aTurb.cval()[i]
										* surfaceNeighbourSide.density[0].cval()[i];

						if (gasPhase.turbulence->bField())
							surfaceNeighbourSide.rhobTurb.val()[i] =
									surfaceNeighbourSide.bTurb.cval()[i]
											* surfaceNeighbourSide.density[0].cval()[i];
					}
				}

				std::valarray<scalar> concentrations_i(
						surfaceNeighbourSide.concentration.v.size());
				for (std::size_t k = 0; k < concentrations_i.size(); ++k)
					concentrations_i[k] =
							surfaceNeighbourSide.concentration.v[k].cval()[i];

				surfaceNeighbourSide.internalEnergy.val()[i] =
						surfaceNeighbourSide.phaseThermodynamics->UvFromp(
								concentrations_i,
								surfaceNeighbourSide.pressure.cval()[i]);

				surfaceNeighbourSide.temperature.val()[i] =
						surfaceNeighbourSide.phaseThermodynamics->TFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy.cval()[i]);

				const scalar v2 { surfaceNeighbourSide.velocity.cval()[i]
						& surfaceNeighbourSide.velocity.cval()[i] };

				surfaceNeighbourSide.totalEnergy.val()[i] =
						surfaceNeighbourSide.internalEnergy.cval()[i]
								+ surfaceNeighbourSide.density[0].cval()[i] * v2
										* 0.5
								+ surfaceNeighbourSide.rhokTurb.cval()[i];

				surfaceNeighbourSide.HelmholtzEnergy.val()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Fv(
								concentrations_i,
								surfaceNeighbourSide.temperature.cval()[i]);

				surfaceNeighbourSide.entropy.val()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Sv(
								concentrations_i,
								surfaceNeighbourSide.temperature.cval()[i]);
			}
	}
	const auto TVDEndTime = std::chrono::high_resolution_clock::now();
	timeForTVD += std::chrono::duration_cast<std::chrono::milliseconds>(
			TVDEndTime - TVDStartTime).count();

	const auto HancockStartTime = std::chrono::high_resolution_clock::now();
	const auto HancockEndTime = std::chrono::high_resolution_clock::now();
	timeForHancock += std::chrono::duration_cast<std::chrono::milliseconds>(
			HancockEndTime - HancockStartTime).count();

	/*Update parallel boundary conditions*/
	parallelism.correctBoundaryValues(surfaceOwnerSide);

	/*Outer neighbour side.*/
	const std::size_t nonExistentCell = mesh_.nonexistCell();
	for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
		if (mesh_.surfaceNeighbour()[i] == nonExistentCell)
		{
			/*Coping from surface owner side.*/
			surfaceNeighbourSide.concentration.v[0].val()[i] = 0;
			surfaceNeighbourSide.density[0].val()[i] = 0;
			for (std::size_t k = 1;
					k < surfaceNeighbourSide.concentration.v.size(); ++k)
			{
				surfaceNeighbourSide.concentration.v[k].val()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.concentration.v[k].cval()[i],
								surfaceOwnerSide.concentration.v[k].boundCond()[i],
								i, i);

				surfaceNeighbourSide.density[k].val()[i] =
						surfaceNeighbourSide.concentration.v[k].cval()[i]
								* surfaceNeighbourSide.phaseThermodynamics->Mv()[k
										- 1];

				surfaceNeighbourSide.concentration.v[0].val()[i] +=
						surfaceNeighbourSide.concentration.v[k].cval()[i];

				surfaceNeighbourSide.density[0].val()[i] +=
						surfaceNeighbourSide.density[k].cval()[i];
			}

			surfaceNeighbourSide.velocity.val()[i] =
					boundaryConditionValueCalc.boundaryConditionValueCell(
							surfaceOwnerSide.velocity.cval()[i],
							surfaceOwnerSide.velocity.boundCond()[i], i, i);

			surfaceNeighbourSide.momentum.val()[i] =
					surfaceNeighbourSide.velocity.cval()[i]
							* surfaceNeighbourSide.density[0].cval()[i];

			surfaceNeighbourSide.pressure.val()[i] =
					boundaryConditionValueCalc.boundaryConditionValueCell(
							surfaceOwnerSide.pressure.cval()[i],
							surfaceOwnerSide.pressure.boundCond()[i], i, i);

			if (surfaceNeighbourSide.turbulence->turbulence())
			{
				surfaceNeighbourSide.kTurb.val()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.kTurb.cval()[i],
								surfaceOwnerSide.kTurb.boundCond()[i], i, i);

				surfaceNeighbourSide.rhokTurb.val()[i] =
						surfaceNeighbourSide.density[0].cval()[i]
								* surfaceNeighbourSide.kTurb.cval()[i];

				surfaceNeighbourSide.epsTurb.val()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.epsTurb.cval()[i],
								surfaceOwnerSide.epsTurb.boundCond()[i], i, i);

				surfaceNeighbourSide.rhoepsTurb.val()[i] =
						surfaceNeighbourSide.density[0].cval()[i]
								* surfaceNeighbourSide.epsTurb.cval()[i];

				if (surfaceNeighbourSide.turbulence->aField())
				{
					surfaceNeighbourSide.aTurb.val()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.aTurb.cval()[i],
									surfaceOwnerSide.aTurb.boundCond()[i], i,
									i);

					surfaceNeighbourSide.rhoaTurb.val()[i] =
							surfaceNeighbourSide.aTurb.cval()[i]
									* surfaceNeighbourSide.density[0].cval()[i];

					if (surfaceNeighbourSide.turbulence->bField())
					{
						surfaceNeighbourSide.bTurb.val()[i] =
								boundaryConditionValueCalc.boundaryConditionValueCell(
										surfaceOwnerSide.bTurb.cval()[i],
										surfaceOwnerSide.bTurb.boundCond()[i],
										i, i);

						surfaceNeighbourSide.rhobTurb.val()[i] =
								surfaceNeighbourSide.density[0].cval()[i]
										* surfaceNeighbourSide.bTurb.cval()[i];
					}
				}
			}

			std::valarray<scalar> concentrations_i(
					surfaceNeighbourSide.concentration.v.size());
			for (std::size_t k = 0; k < concentrations_i.size(); ++k)
				concentrations_i[k] =
						surfaceNeighbourSide.concentration.v[k].cval()[i];

			surfaceNeighbourSide.internalEnergy.val()[i] =
					surfaceNeighbourSide.phaseThermodynamics->UvFromp(
							concentrations_i,
							surfaceNeighbourSide.pressure.cval()[i]);

			surfaceNeighbourSide.temperature.val()[i] =
					surfaceNeighbourSide.phaseThermodynamics->TFromUv(
							concentrations_i,
							surfaceNeighbourSide.internalEnergy.cval()[i]);

			const scalar v2 { surfaceNeighbourSide.velocity.cval()[i]
					& surfaceNeighbourSide.velocity.cval()[i] };

			surfaceNeighbourSide.totalEnergy.val()[i] =
					surfaceNeighbourSide.internalEnergy.cval()[i]
							+ surfaceNeighbourSide.density[0].cval()[i] * v2
									* 0.5
							+ surfaceNeighbourSide.rhokTurb.cval()[i];

			surfaceNeighbourSide.HelmholtzEnergy.val()[i] =
					surfaceNeighbourSide.phaseThermodynamics->Fv(
							concentrations_i,
							surfaceNeighbourSide.temperature.cval()[i]);

			surfaceNeighbourSide.entropy.val()[i] =
					surfaceNeighbourSide.phaseThermodynamics->Sv(
							concentrations_i,
							surfaceNeighbourSide.temperature.cval()[i]);
		}

	/*Numerical fluxes.*/
	const auto FlowCalculationStartTime =
			std::chrono::high_resolution_clock::now();

	/*Creating fields of flux solver.*/
	const auto [NumFluxFlows, star] = fsolver.calculateFlows(surfaceOwnerSide,
			surfaceNeighbourSide);

	const auto FlowCalculationEndTime =
			std::chrono::high_resolution_clock::now();
	timeForFlowCalculation += std::chrono::duration_cast<
			std::chrono::milliseconds>(
			FlowCalculationEndTime - FlowCalculationStartTime).count();

	/*Time integration.*/
	const auto TimeIntegrationStartTime =
			std::chrono::high_resolution_clock::now();
	{
		const scalar timestep = mesh_.timestep();

		gasPhase.totalEnergy -= divergence(NumFluxFlows.totalEnergy) * timestep;

		gasPhase.momentum -= divergence(NumFluxFlows.momentum) * timestep;

		if (gravitation.first)
		{
			gasPhase.totalEnergy += gasPhase.density[0]
					* (gasPhase.velocity & gravitation.second) * timestep;

			gasPhase.momentum += gasPhase.density[0] * gravitation.second
					* timestep;
		}

		for (std::size_t k = 0; k < gasPhase.density.size(); ++k)
			gasPhase.density[k] -= divergence(NumFluxFlows.density[k])
					* timestep;

		if (gasPhase.turbulence->turbulence())
		{
			gasPhase.rhokTurb -= divergence(NumFluxFlows.rhokTurb) * timestep;

			gasPhase.rhoepsTurb -= divergence(NumFluxFlows.rhoepsTurb)
					* timestep;

			if (gasPhase.turbulence->aField())
			{
				gasPhase.rhoaTurb -= divergence(NumFluxFlows.rhoaTurb)
						* timestep;

				if (gasPhase.turbulence->bField())
					gasPhase.rhobTurb -= divergence(NumFluxFlows.rhobTurb)
							* timestep;
			}
		}

		gasPhase.concentration.v[0].val() = 0;
		for (std::size_t k = 1; k < gasPhase.concentration.v.size(); ++k)
		{
			gasPhase.concentration.v[k].val() = (gasPhase.density[k]
					/ gasPhase.phaseThermodynamics->Mv()[k - 1]).cval();

			gasPhase.concentration.v[0] += gasPhase.concentration.v[k];
		}

		gasPhase.velocity.val() =
				(gasPhase.momentum / gasPhase.density[0]).cval();

		{
			const auto v2 = gasPhase.velocity & gasPhase.velocity;

			gasPhase.internalEnergy.val() =
					(gasPhase.totalEnergy - gasPhase.density[0] * v2 * 0.5
							- gasPhase.rhokTurb).cval();
		}

		gasPhase.pressure.val() = gasPhase.phaseThermodynamics->pFromUv(
				gasPhase.concentration.p, gasPhase.internalEnergy.cval());

		gasPhase.temperature.val() = gasPhase.phaseThermodynamics->TFromUv(
				gasPhase.concentration.p, gasPhase.internalEnergy.cval());

		gasPhase.HelmholtzEnergy.val() = gasPhase.phaseThermodynamics->Fv(
				gasPhase.concentration.p, gasPhase.temperature.cval());

		gasPhase.entropy.val() = gasPhase.phaseThermodynamics->Sv(
				gasPhase.concentration.p, gasPhase.temperature.cval());

		if (gasPhase.turbulence->turbulence())
		{
			gasPhase.kTurb.val() =
					(gasPhase.rhokTurb / gasPhase.density[0]).cval();

			gasPhase.epsTurb.val() =
					(gasPhase.rhoepsTurb / gasPhase.density[0]).cval();

			if (gasPhase.turbulence->aField())
			{
				gasPhase.aTurb.val() =
						(gasPhase.rhoaTurb / gasPhase.density[0]).cval();

				if (gasPhase.turbulence->bField())
					gasPhase.bTurb.val() = (gasPhase.rhobTurb
							/ gasPhase.density[0]).cval();
			}
		}
	}
	const auto TimeIntegrationEndTime =
			std::chrono::high_resolution_clock::now();
	timeForTimeIntegration += std::chrono::duration_cast<
			std::chrono::milliseconds>(
			TimeIntegrationEndTime - TimeIntegrationStartTime).count();

	return star;
}
