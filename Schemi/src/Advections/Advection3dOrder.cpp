/*
 * Advection3dOrder.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "Advection3dOrder.hpp"

#include <chrono>

#include "divergence.hpp"
#include "gradient.hpp"
#include "thirdOrderAveraging.hpp"

namespace schemi
{
class abstractLimiter;
} /* namespace schemi */

schemi::starFields schemi::Advection3dOrder(
		homogeneousPhase<cubicCell> & gasPhase, const abstractLimiter & limiter,
		const abstractFlowSolver & fsolver, std::pair<bool, vector> gravitation,
		const boundaryConditionValue & boundaryConditionValueCalc,
		scalar & timeForTVD, scalar & timeForHancock,
		scalar & timeForFlowCalculation, scalar & timeForTimeIntegration,
		const MPIHandler & parallelism)
{
	auto & mesh { gasPhase.pressure.meshRef() };

	/*Creating surface fields.*/
	homogeneousPhase<quadraticSurface> surfaceOwnerSide { bunchOfFields<
			quadraticSurface>(gasPhase),
			transportCoefficients<quadraticSurface>(mesh),
			gasPhase.phaseThermodynamics, gasPhase.turbulenceSources,
			gasPhase.transportModel };

	homogeneousPhase<quadraticSurface> surfaceNeighbourSide { bunchOfFields<
			quadraticSurface>(gasPhase),
			transportCoefficients<quadraticSurface>(mesh),
			gasPhase.phaseThermodynamics, gasPhase.turbulenceSources,
			gasPhase.transportModel };

	/*Creating limited gradients.*/
	std::vector<volumeField<std::vector<vector>> > concentrationTVDGradient {
			gasPhase.phaseThermodynamics->Mv().size(), volumeField<
					std::vector<vector>>(mesh, std::vector<vector>(0)) };
	volumeField<std::vector<tensor>> velocityTVDGradient { mesh, std::vector<
			tensor>(0) };
	volumeField<std::vector<vector>> pressureTVDGradient { mesh, std::vector<
			vector>(0) };
	volumeField<std::vector<vector>> kTVDGradient { mesh, std::vector<vector>(0) };
	volumeField<std::vector<vector>> epsilonTVDGradient { mesh, std::vector<
			vector>(0) };
	volumeField<std::vector<tensor>> aTVDGradient { mesh, std::vector<tensor>(0) };
	volumeField<std::vector<vector>> bTVDGradient { mesh, std::vector<vector>(0) };

	/*TVD Reconstruction.*/
	const auto TVDStartTime = std::chrono::high_resolution_clock::now();
	{
		/*TVD limiters.*/
		for (std::size_t k = 0; k < concentrationTVDGradient.size(); ++k)
		{
			const auto concentrationGradient_k = grad(
					gasPhase.concentration.v[k + 1], boundaryConditionValueCalc,
					k + 1);

			const auto surfaceConcentarionGradient_k =
					gradientLinearInterpolate(concentrationGradient_k,
							gasPhase.concentration.v[k + 1],
							boundaryConditionValueCalc, k + 1);

			concentrationTVDGradient[k] = thirdOrderAveraging(
					surfaceConcentarionGradient_k, limiter);
		}

		const auto velocityGradient = grad(gasPhase.velocity,
				boundaryConditionValueCalc);

		const auto surfaceVelocityGradient = gradientLinearInterpolate(
				velocityGradient, gasPhase.velocity,
				boundaryConditionValueCalc);

		velocityTVDGradient = thirdOrderAveraging(surfaceVelocityGradient,
				limiter);

		const auto pressureGradient = grad(gasPhase.pressure,
				boundaryConditionValueCalc);

		const auto surfacePressureGradient = gradientLinearInterpolate(
				pressureGradient, gasPhase.pressure,
				boundaryConditionValueCalc);

		pressureTVDGradient = thirdOrderAveraging(surfacePressureGradient,
				limiter);

		if (gasPhase.turbulenceSources->turbulence)
		{
			const auto kGradient = grad(gasPhase.kTurb,
					boundaryConditionValueCalc);

			const auto surfacekGradient = gradientLinearInterpolate(kGradient,
					gasPhase.kTurb, boundaryConditionValueCalc);

			kTVDGradient = thirdOrderAveraging(surfacekGradient, limiter);

			const auto epsilonGradient = grad(gasPhase.epsTurb,
					boundaryConditionValueCalc);

			const auto surfaceepsilonGradient = gradientLinearInterpolate(
					epsilonGradient, gasPhase.epsTurb,
					boundaryConditionValueCalc);

			epsilonTVDGradient = thirdOrderAveraging(surfaceepsilonGradient,
					limiter);

			if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::BHRKLSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::kEpsASource))
			{
				const auto aGradient = grad(gasPhase.aTurb,
						boundaryConditionValueCalc);

				const auto surfaceaGradient = gradientLinearInterpolate(
						aGradient, gasPhase.aTurb, boundaryConditionValueCalc);

				aTVDGradient = thirdOrderAveraging(surfaceaGradient, limiter);

				if ((gasPhase.turbulenceSources->model
						== turbulenceModel::BHRSource)
						|| (gasPhase.turbulenceSources->model
								== turbulenceModel::BHRKLSource))
				{
					const auto bGradient = grad(gasPhase.bTurb,
							boundaryConditionValueCalc);

					const auto surfacebGradient = gradientLinearInterpolate(
							bGradient, gasPhase.bTurb,
							boundaryConditionValueCalc);

					bTVDGradient = thirdOrderAveraging(surfacebGradient,
							limiter);
				}
			}
		}

		for (std::size_t i = 0; i < mesh.cellsSize(); ++i)
		{
			const std::vector<std::size_t> & surfacesOfCell_i {
					mesh.surfacesOfCells()[i] };

			for (std::size_t j = 0; j < surfacesOfCell_i.size(); ++j)
			{
				const std::size_t surfaceIndex { surfacesOfCell_i[j] };

				const vector deltaVector = mesh.surfaces()[surfaceIndex].rC()
						- mesh.cells()[i].rC();

				std::valarray<scalar> reconstructedConcentrationsValues(
						gasPhase.phaseThermodynamics->Mv().size());

				for (std::size_t k = 0;
						k < reconstructedConcentrationsValues.size(); ++k)
					reconstructedConcentrationsValues[k] =
							(concentrationTVDGradient[k].ref()[i][j]
									& deltaVector)
									+ gasPhase.concentration.v[k + 1].ref()[i];

				const vector reconstructedVelocityValue =
						(velocityTVDGradient.ref()[i][j] & deltaVector)
								+ gasPhase.velocity.ref()[i];

				const scalar reconstructedPressureValue {
						(pressureTVDGradient.ref()[i][j] & deltaVector)
								+ gasPhase.pressure.ref()[i] };

				scalar reconstructedkValue { 0 };
				scalar reconstructedEpsilon { 0 };
				vector reconstructedaValue { 0 };
				scalar reconstructedbValue { 0 };
				if (gasPhase.turbulenceSources->turbulence)
				{
					reconstructedkValue = (kTVDGradient.ref()[i][j]
							& deltaVector) + gasPhase.kTurb.ref()[i];

					reconstructedEpsilon = (epsilonTVDGradient.ref()[i][j]
							& deltaVector) + gasPhase.epsTurb.ref()[i];

					if ((gasPhase.turbulenceSources->model
							== turbulenceModel::BHRSource)
							|| (gasPhase.turbulenceSources->model
									== turbulenceModel::BHRKLSource)
							|| (gasPhase.turbulenceSources->model
									== turbulenceModel::kEpsASource))
					{
						reconstructedaValue = (aTVDGradient.ref()[i][j]
								& deltaVector) + gasPhase.aTurb.ref()[i];

						if ((gasPhase.turbulenceSources->model
								== turbulenceModel::BHRSource)
								|| (gasPhase.turbulenceSources->model
										== turbulenceModel::BHRKLSource))
							reconstructedbValue = (bTVDGradient.ref()[i][j]
									& deltaVector) + gasPhase.bTurb.ref()[i];
					}
				}

				if (mesh.surfaceOwner()[surfaceIndex] == i)
				{
					for (std::size_t k = 0;
							k < reconstructedConcentrationsValues.size(); ++k)
						surfaceOwnerSide.concentration.v[k + 1].ref_r()[surfaceIndex] =
								reconstructedConcentrationsValues[k];

					surfaceOwnerSide.velocity.ref_r()[surfaceIndex] =
							reconstructedVelocityValue;

					surfaceOwnerSide.pressure.ref_r()[surfaceIndex] =
							reconstructedPressureValue;

					if (gasPhase.turbulenceSources->turbulence)
					{
						surfaceOwnerSide.kTurb.ref_r()[surfaceIndex] =
								reconstructedkValue;

						surfaceOwnerSide.epsTurb.ref_r()[surfaceIndex] =
								reconstructedEpsilon;

						if ((gasPhase.turbulenceSources->model
								== turbulenceModel::BHRSource)
								|| (gasPhase.turbulenceSources->model
										== turbulenceModel::BHRKLSource)
								|| (gasPhase.turbulenceSources->model
										== turbulenceModel::kEpsASource))
						{
							surfaceOwnerSide.aTurb.ref_r()[surfaceIndex] =
									reconstructedaValue;

							if ((gasPhase.turbulenceSources->model
									== turbulenceModel::BHRSource)
									|| (gasPhase.turbulenceSources->model
											== turbulenceModel::BHRKLSource))
								surfaceOwnerSide.bTurb.ref_r()[surfaceIndex] =
										reconstructedbValue;
						}
					}
				}
				else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				{
					for (std::size_t k = 0;
							k < reconstructedConcentrationsValues.size(); ++k)
						surfaceNeighbourSide.concentration.v[k + 1].ref_r()[surfaceIndex] =
								reconstructedConcentrationsValues[k];

					surfaceNeighbourSide.velocity.ref_r()[surfaceIndex] =
							reconstructedVelocityValue;

					surfaceNeighbourSide.pressure.ref_r()[surfaceIndex] =
							reconstructedPressureValue;

					if (gasPhase.turbulenceSources->turbulence)
					{
						surfaceNeighbourSide.kTurb.ref_r()[surfaceIndex] =
								reconstructedkValue;

						surfaceNeighbourSide.epsTurb.ref_r()[surfaceIndex] =
								reconstructedEpsilon;

						if ((gasPhase.turbulenceSources->model
								== turbulenceModel::BHRSource)
								|| (gasPhase.turbulenceSources->model
										== turbulenceModel::BHRKLSource)
								|| (gasPhase.turbulenceSources->model
										== turbulenceModel::kEpsASource))
						{
							surfaceNeighbourSide.aTurb.ref_r()[surfaceIndex] =
									reconstructedaValue;

							if ((gasPhase.turbulenceSources->model
									== turbulenceModel::BHRSource)
									|| (gasPhase.turbulenceSources->model
											== turbulenceModel::BHRKLSource))
								surfaceNeighbourSide.bTurb.ref_r()[surfaceIndex] =
										reconstructedbValue;
						}
					}
				}
				else
					throw exception("Couldn't choose side to add",
							errors::systemError);
			}
		}

		/*Recalculation of other quantities for owner side*/
		surfaceOwnerSide.concentration.v[0].ref_r() = 0;
		surfaceOwnerSide.density[0].ref_r() = 0;
		for (std::size_t k = 1; k < surfaceOwnerSide.concentration.v.size();
				++k)
		{
			surfaceOwnerSide.concentration.v[0].ref_r() +=
					surfaceOwnerSide.concentration.v[k].ref();

			surfaceOwnerSide.density[k].ref_r() =
					surfaceOwnerSide.concentration.v[k].ref()
							* surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			surfaceOwnerSide.density[0].ref_r() +=
					surfaceOwnerSide.density[k].ref();
		}

		if (gasPhase.turbulenceSources->turbulence)
		{
			surfaceOwnerSide.rhokTurb.ref_r() = astProduct(
					surfaceOwnerSide.kTurb, surfaceOwnerSide.density[0]).ref();

			surfaceOwnerSide.rhoepsTurb.ref_r() = surfaceOwnerSide.epsTurb.ref()
					* surfaceOwnerSide.density[0].ref();

			if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::BHRKLSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::kEpsASource))
			{
				surfaceOwnerSide.rhoaTurb.ref_r() =
						astProduct(surfaceOwnerSide.aTurb,
								surfaceOwnerSide.density[0]).ref();

				if ((gasPhase.turbulenceSources->model
						== turbulenceModel::BHRSource)
						|| (gasPhase.turbulenceSources->model
								== turbulenceModel::BHRKLSource))
					surfaceOwnerSide.rhobTurb.ref_r() =
							surfaceOwnerSide.bTurb.ref()
									* surfaceOwnerSide.density[0].ref();
			}
		}

		surfaceOwnerSide.internalEnergy.ref_r() =
				surfaceOwnerSide.phaseThermodynamics->UvFromp(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.pressure.ref());

		surfaceOwnerSide.temperature.ref_r() =
				surfaceOwnerSide.phaseThermodynamics->TFromUv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.internalEnergy.ref());

		surfaceOwnerSide.momentum.ref_r() = astProduct(
				surfaceOwnerSide.velocity, surfaceOwnerSide.density[0]).ref();

		{
			const auto v2 = ampProduct(surfaceOwnerSide.velocity,
					surfaceOwnerSide.velocity);

			surfaceOwnerSide.totalEnergy.ref_r() =
					surfaceOwnerSide.internalEnergy.ref()
							+ surfaceOwnerSide.density[0].ref() * v2.ref() * 0.5
							+ surfaceOwnerSide.rhokTurb.ref();
		}

		surfaceOwnerSide.HelmholtzEnergy.ref_r() =
				surfaceOwnerSide.phaseThermodynamics->Fv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature.ref());

		surfaceOwnerSide.entropy.ref_r() =
				surfaceOwnerSide.phaseThermodynamics->Sv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature.ref());

		/*Recalculation of other quantities for neighbour side*/
		const std::size_t nonExistentCell = mesh.nonexistCell();
		for (std::size_t i = 0; i < mesh.surfacesSize(); ++i)
			if (mesh.surfaceNeighbour()[i] != nonExistentCell)
			{
				surfaceNeighbourSide.concentration.v[0].ref_r()[i] = 0;
				surfaceNeighbourSide.density[0].ref_r()[i] = 0;
				for (std::size_t k = 1;
						k < surfaceNeighbourSide.concentration.v.size(); ++k)
				{
					surfaceNeighbourSide.concentration.v[0].ref_r()[i] +=
							surfaceNeighbourSide.concentration.v[k].ref()[i];

					surfaceNeighbourSide.density[k].ref_r()[i] =
							surfaceNeighbourSide.concentration.v[k].ref()[i]
									* gasPhase.phaseThermodynamics->Mv()[k - 1];

					surfaceNeighbourSide.density[0].ref_r()[i] +=
							surfaceNeighbourSide.density[k].ref()[i];
				}

				surfaceNeighbourSide.momentum.ref_r()[i] =
						surfaceNeighbourSide.velocity.ref()[i]
								* surfaceNeighbourSide.density[0].ref()[i];

				if (gasPhase.turbulenceSources->turbulence)
				{
					surfaceNeighbourSide.rhokTurb.ref_r()[i] =
							surfaceNeighbourSide.kTurb.ref()[i]
									* surfaceNeighbourSide.density[0].ref()[i];

					surfaceNeighbourSide.rhoepsTurb.ref_r()[i] =
							surfaceNeighbourSide.epsTurb.ref()[i]
									* surfaceNeighbourSide.density[0].ref()[i];

					if ((gasPhase.turbulenceSources->model
							== turbulenceModel::BHRSource)
							|| (gasPhase.turbulenceSources->model
									== turbulenceModel::BHRKLSource)
							|| (gasPhase.turbulenceSources->model
									== turbulenceModel::kEpsASource))
					{
						surfaceNeighbourSide.rhoaTurb.ref_r()[i] =
								surfaceNeighbourSide.aTurb.ref()[i]
										* surfaceNeighbourSide.density[0].ref()[i];

						if ((gasPhase.turbulenceSources->model
								== turbulenceModel::BHRSource)
								|| (gasPhase.turbulenceSources->model
										== turbulenceModel::BHRKLSource))
							surfaceNeighbourSide.rhobTurb.ref_r()[i] =
									surfaceNeighbourSide.bTurb.ref()[i]
											* surfaceNeighbourSide.density[0].ref()[i];
					}
				}

				std::valarray<scalar> concentrations_i(
						surfaceNeighbourSide.concentration.v.size());
				for (std::size_t k = 0; k < concentrations_i.size(); ++k)
					concentrations_i[k] =
							surfaceNeighbourSide.concentration.v[k].ref()[i];

				surfaceNeighbourSide.internalEnergy.ref_r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->UvFromp(
								concentrations_i,
								surfaceNeighbourSide.pressure.ref()[i]);

				surfaceNeighbourSide.temperature.ref_r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->TFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy.ref()[i]);

				const scalar v2 { surfaceNeighbourSide.velocity.ref()[i]
						& surfaceNeighbourSide.velocity.ref()[i] };

				surfaceNeighbourSide.totalEnergy.ref_r()[i] =
						surfaceNeighbourSide.internalEnergy.ref()[i]
								+ surfaceNeighbourSide.density[0].ref()[i] * v2
										* 0.5
								+ surfaceNeighbourSide.rhokTurb.ref()[i];

				surfaceNeighbourSide.HelmholtzEnergy.ref_r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Fv(
								concentrations_i,
								surfaceNeighbourSide.temperature.ref()[i]);

				surfaceNeighbourSide.entropy.ref_r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Sv(
								concentrations_i,
								surfaceNeighbourSide.temperature.ref()[i]);
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
	const std::size_t nonExistentCell = mesh.nonexistCell();
	for (std::size_t i = 0; i < mesh.surfacesSize(); ++i)
		if (mesh.surfaceNeighbour()[i] == nonExistentCell)
		{
			/*Coping from surface owner side.*/
			surfaceNeighbourSide.concentration.v[0].ref_r()[i] = 0;
			surfaceNeighbourSide.density[0].ref_r()[i] = 0;
			for (std::size_t k = 1;
					k < surfaceNeighbourSide.concentration.v.size(); ++k)
			{
				surfaceNeighbourSide.concentration.v[k].ref_r()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.concentration.v[k].ref()[i],
								surfaceOwnerSide.concentration.v[k].boundCond()[i],
								i, i, k);

				surfaceNeighbourSide.density[k].ref_r()[i] =
						surfaceNeighbourSide.concentration.v[k].ref()[i]
								* surfaceNeighbourSide.phaseThermodynamics->Mv()[k
										- 1];

				surfaceNeighbourSide.concentration.v[0].ref_r()[i] +=
						surfaceNeighbourSide.concentration.v[k].ref()[i];

				surfaceNeighbourSide.density[0].ref_r()[i] +=
						surfaceNeighbourSide.density[k].ref()[i];
			}

			surfaceNeighbourSide.velocity.ref_r()[i] =
					boundaryConditionValueCalc.boundaryConditionValueCell(
							surfaceOwnerSide.velocity.ref()[i],
							surfaceOwnerSide.velocity.boundCond()[i], i, i);

			surfaceNeighbourSide.momentum.ref_r()[i] =
					surfaceNeighbourSide.velocity.ref_r()[i]
							* surfaceNeighbourSide.density[0].ref()[i];

			surfaceNeighbourSide.pressure.ref_r()[i] =
					boundaryConditionValueCalc.boundaryConditionValueCell(
							surfaceOwnerSide.pressure.ref()[i],
							surfaceOwnerSide.pressure.boundCond()[i], i, i);

			if (surfaceNeighbourSide.turbulenceSources->turbulence)
			{
				surfaceNeighbourSide.kTurb.ref_r()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.kTurb.ref()[i],
								surfaceOwnerSide.kTurb.boundCond()[i], i, i);

				surfaceNeighbourSide.rhokTurb.ref_r()[i] =
						surfaceNeighbourSide.density[0].ref_r()[i]
								* surfaceNeighbourSide.kTurb.ref()[i];

				surfaceNeighbourSide.epsTurb.ref_r()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.epsTurb.ref()[i],
								surfaceOwnerSide.epsTurb.boundCond()[i], i, i);

				surfaceNeighbourSide.rhoepsTurb.ref_r()[i] =
						surfaceNeighbourSide.density[0].ref()[i]
								* surfaceNeighbourSide.epsTurb.ref()[i];

				if ((surfaceNeighbourSide.turbulenceSources->model
						== turbulenceModel::BHRSource)
						|| (surfaceNeighbourSide.turbulenceSources->model
								== turbulenceModel::BHRKLSource)
						|| (surfaceNeighbourSide.turbulenceSources->model
								== turbulenceModel::kEpsASource))
				{
					surfaceNeighbourSide.aTurb.ref_r()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.aTurb.ref()[i],
									surfaceOwnerSide.aTurb.boundCond()[i], i,
									i);

					surfaceNeighbourSide.rhoaTurb.ref_r()[i] =
							surfaceNeighbourSide.aTurb.ref_r()[i]
									* surfaceNeighbourSide.density[0].ref()[i];

					if ((surfaceNeighbourSide.turbulenceSources->model
							== turbulenceModel::BHRSource)
							|| (surfaceNeighbourSide.turbulenceSources->model
									== turbulenceModel::BHRKLSource))
					{
						surfaceNeighbourSide.bTurb.ref_r()[i] =
								boundaryConditionValueCalc.boundaryConditionValueCell(
										surfaceOwnerSide.bTurb.ref()[i],
										surfaceOwnerSide.bTurb.boundCond()[i],
										i, i);

						surfaceNeighbourSide.rhobTurb.ref_r()[i] =
								surfaceNeighbourSide.density[0].ref_r()[i]
										* surfaceNeighbourSide.bTurb.ref()[i];
					}
				}
			}

			std::valarray<scalar> concentrations_i(
					surfaceNeighbourSide.concentration.v.size());
			for (std::size_t k = 0; k < concentrations_i.size(); ++k)
				concentrations_i[k] =
						surfaceNeighbourSide.concentration.v[k].ref()[i];

			surfaceNeighbourSide.internalEnergy.ref_r()[i] =
					surfaceNeighbourSide.phaseThermodynamics->UvFromp(
							concentrations_i,
							surfaceNeighbourSide.pressure.ref()[i]);

			surfaceNeighbourSide.temperature.ref_r()[i] =
					surfaceNeighbourSide.phaseThermodynamics->TFromUv(
							concentrations_i,
							surfaceNeighbourSide.internalEnergy.ref()[i]);

			const scalar v2 { surfaceNeighbourSide.velocity.ref()[i]
					& surfaceNeighbourSide.velocity.ref()[i] };

			surfaceNeighbourSide.totalEnergy.ref_r()[i] =
					surfaceNeighbourSide.internalEnergy.ref_r()[i]
							+ surfaceNeighbourSide.density[0].ref()[i] * v2
									* 0.5
							+ surfaceNeighbourSide.rhokTurb.ref()[i];

			surfaceNeighbourSide.HelmholtzEnergy.ref_r()[i] =
					surfaceNeighbourSide.phaseThermodynamics->Fv(
							concentrations_i,
							surfaceNeighbourSide.temperature.ref()[i]);

			surfaceNeighbourSide.entropy.ref_r()[i] =
					surfaceNeighbourSide.phaseThermodynamics->Sv(
							concentrations_i,
							surfaceNeighbourSide.temperature.ref()[i]);
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
		const scalar timestep = mesh.timestep();

		gasPhase.totalEnergy.ref_r() -=
				divergence(NumFluxFlows.totalEnergy).ref() * timestep;

		gasPhase.momentum.ref_r() -= astProduct(
				divergence(NumFluxFlows.momentum), timestep).ref();

		if (gravitation.first)
		{
			gasPhase.totalEnergy.ref_r() += gasPhase.density[0].ref()
					* ampProduct(gasPhase.velocity, gravitation.second).ref()
					* timestep;

			gasPhase.momentum.ref_r() += astProduct(
					astProduct(gasPhase.density[0], gravitation.second),
					timestep).ref();
		}

		for (std::size_t k = 0; k < gasPhase.density.size(); ++k)
		{
			gasPhase.density[k].ref_r() -=
					divergence(NumFluxFlows.density[k]).ref() * timestep;

			if (k != 0)
				std::replace_if(std::begin(gasPhase.density[k].ref_r()),
						std::end(gasPhase.density[k].ref_r()),
						[](const scalar value) 
						{
							return value < 0.0;
						}, 0.0);
		}

		for (std::size_t i = 0; i < gasPhase.density[0].size(); ++i)
		{
			scalar sumDensity { 0. };

			for (std::size_t k = 1; k < gasPhase.density.size(); ++k)
				sumDensity += gasPhase.density[k].ref()[i];

			if (sumDensity > gasPhase.density[0].ref()[i])
				gasPhase.density[0].ref_r()[i] = sumDensity;
		}

		if (gasPhase.turbulenceSources->turbulence)
		{
			gasPhase.rhokTurb.ref_r() -= divergence(NumFluxFlows.rhokTurb).ref()
					* timestep;

			gasPhase.rhoepsTurb.ref_r() -=
					divergence(NumFluxFlows.rhoepsTurb).ref() * timestep;

			std::replace_if(std::begin(gasPhase.rhokTurb.ref_r()),
					std::end(gasPhase.rhokTurb.ref_r()),
					[&gasPhase](
							const scalar value) 
							{
								return value < gasPhase.turbulenceSources->turbPar->mink();
							}, gasPhase.turbulenceSources->turbPar->mink());

			std::replace_if(std::begin(gasPhase.rhoepsTurb.ref_r()),
					std::end(gasPhase.rhoepsTurb.ref_r()),
					[&gasPhase](
							const scalar value) 
							{
								return value < gasPhase.turbulenceSources->turbPar->mineps();
							}, gasPhase.turbulenceSources->turbPar->mineps());

			if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::BHRKLSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::kEpsASource))
			{
				gasPhase.rhoaTurb.ref_r() -= astProduct(
						divergence(NumFluxFlows.rhoaTurb), timestep).ref();

				if ((gasPhase.turbulenceSources->model
						== turbulenceModel::BHRSource)
						|| (gasPhase.turbulenceSources->model
								== turbulenceModel::BHRKLSource))
				{
					gasPhase.rhobTurb.ref_r() -= divergence(
							NumFluxFlows.rhobTurb).ref() * timestep;

					std::replace_if(std::begin(gasPhase.rhobTurb.ref_r()),
							std::end(gasPhase.rhobTurb.ref_r()),
							[&gasPhase](
									const scalar value) 
									{
										return value < gasPhase.turbulenceSources->turbPar->minb_value;
									},
							gasPhase.turbulenceSources->turbPar->minb_value);
				}
			}
		}

		gasPhase.concentration.v[0].ref_r() = 0;
		for (std::size_t k = 1; k < gasPhase.concentration.v.size(); ++k)
		{
			gasPhase.concentration.v[k].ref_r() = gasPhase.density[k].ref()
					/ gasPhase.phaseThermodynamics->Mv()[k - 1];

			gasPhase.concentration.v[0].ref_r() +=
					gasPhase.concentration.v[k].ref();
		}

		gasPhase.velocity.ref_r() = division(gasPhase.momentum,
				gasPhase.density[0]).ref();

		{
			const auto v2 = ampProduct(gasPhase.velocity, gasPhase.velocity);

			gasPhase.internalEnergy.ref_r() = gasPhase.totalEnergy.ref()
					- gasPhase.density[0].ref() * v2.ref() * 0.5
					- gasPhase.rhokTurb.ref();
		}

		gasPhase.pressure.ref_r() = gasPhase.phaseThermodynamics->pFromUv(
				gasPhase.concentration.p, gasPhase.internalEnergy.ref());

		gasPhase.temperature.ref_r() = gasPhase.phaseThermodynamics->TFromUv(
				gasPhase.concentration.p, gasPhase.internalEnergy.ref());

		gasPhase.HelmholtzEnergy.ref_r() = gasPhase.phaseThermodynamics->Fv(
				gasPhase.concentration.p, gasPhase.temperature.ref());

		gasPhase.entropy.ref_r() = gasPhase.phaseThermodynamics->Sv(
				gasPhase.concentration.p, gasPhase.temperature.ref());

		if (gasPhase.turbulenceSources->turbulence)
		{
			gasPhase.kTurb.ref_r() = gasPhase.rhokTurb.ref()
					/ gasPhase.density[0].ref();

			gasPhase.epsTurb.ref_r() = gasPhase.rhoepsTurb.ref()
					/ gasPhase.density[0].ref();

			if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::BHRKLSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::kEpsASource))
			{
				gasPhase.aTurb.ref_r() = division(gasPhase.rhoaTurb,
						gasPhase.density[0]).ref();

				if ((gasPhase.turbulenceSources->model
						== turbulenceModel::BHRSource)
						|| (gasPhase.turbulenceSources->model
								== turbulenceModel::BHRKLSource))
					gasPhase.bTurb.ref_r() = gasPhase.rhobTurb.ref()
							/ gasPhase.density[0].ref();
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
