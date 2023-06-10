/*
 * Diffusion.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "Diffusion.hpp"

#include <chrono>

#include "harmonicInterpolateScalar.hpp"
#include "divergence.hpp"
#include "fieldProducts.hpp"
#include "gradient.hpp"
#include "linearInterpolate.hpp"

void schemi::Diffusion(homogeneousPhase<cubicCell> & gasPhase,
		const abstractMatrixSolver & msolver,
		const abstractMatrixSolver & msolverForEnthalpy,
		const std::pair<scalar, scalar> & timestepCoeffs,
		scalar & timeForDiffusion,
		const std::vector<boundaryConditionType> & commonConditions,
		const starFields & star, const enthalpyFlowEnum enthalpySolverFlag,
		const bool linearRec, const boundaryConditionValue & bncCalc,
		const volumeField<scalar> & minimalLengthScale,
		const MPIHandler & parallelism, const timestepEnum sourceTimeFlag,
		const bool molMassDiffusionFlag)
{
	const auto DiffusionStartTime { std::chrono::high_resolution_clock::now() };

	auto & mesh_ { gasPhase.pressure.meshRef() };

	const scalar timestep = mesh_.timestep();

	/*Fields for diffusion stage.*/
	surfaceField<vector> entExplDiffFlow { mesh_, vector(0) };
	volumeField<scalar> intermediateTemperature(gasPhase.temperature);
	parallelism.correctBoundaryConditions(intermediateTemperature);
	surfaceField<scalar> surfaceRho { mesh_, 0 };
	volumeField<scalar> CVOld { mesh_, 0 };
	{
		std::vector<boundaryConditionType> bndCon(commonConditions);
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedCv);
		CVOld = volumeField<scalar>(mesh_, 0, bndCon[0], 0, bndCon[1], 0,
				bndCon[2], 0, bndCon[3], 0, bndCon[4], 0, bndCon[5], 0);
	}
	volumeField<scalar> CCVOld { mesh_, 0 };
	volumeField<scalar> CCVNew { mesh_, 0 };
	volumeField<scalar> nonIdealCorrectionOld { mesh_, 0 };
	volumeField<scalar> nonIdealCorrectionNew { mesh_, 0 };
	volumeField<scalar> NonIdRho { mesh_, 0 };
	{
		std::vector<boundaryConditionType> bndCon(commonConditions);
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedNonidealityCorrectionPerDensity);
		NonIdRho = volumeField<scalar>(mesh_, 0, bndCon[0], 0, bndCon[1], 0,
				bndCon[2], 0, bndCon[3], 0, bndCon[4], 0, bndCon[5], 0);
	}

	volumeField<scalar> CvM { mesh_, 0 };
	{
		std::vector<boundaryConditionType> bndCon(commonConditions);
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedCvM);
		CvM = volumeField<scalar>(mesh_, 0, bndCon[0], 0, bndCon[1], 0,
				bndCon[2], 0, bndCon[3], 0, bndCon[4], 0, bndCon[5], 0);
	}

	volumeField<scalar> oldField { mesh_, 0 };

	diffusiveFields diffFieldsOld { mesh_, gasPhase, commonConditions,
			gasPhase.turbulenceSources->turbulence,
			gasPhase.turbulenceSources->model }, diffFieldsNew(diffFieldsOld);

	effectiveTransportCoefficients<quadraticSurface> effectiveCoeffs { mesh_,
			gasPhase.phaseThermodynamics->Mv().size() };

	/*Creating empty objects for SLE matrix.*/
	std::vector<SLEMatrix> massFractionMatrix {
			diffFieldsOld.massFraction.size(), SLEMatrix(
					std::string("Mass fraction")) };
	SLEMatrix velocityMatrix { std::string("Velocity") };
	SLEMatrix temperatureMatrix { std::string("Temperature") };
	SLEMatrix temperatureEnthMatrix { std::string("Temperature-enthalpy") };
	SLEMatrix kMatrix { std::string("k") };
	SLEMatrix epsMatrix { std::string("epsilon") };
	SLEMatrix aMatrix { std::string("a") };
	SLEMatrix bMatrix { std::string("b") };
	/**/

	/*Turbulent coefficients recalculation.*/
	if (gasPhase.turbulenceSources->turbulence)
		gasPhase.recalculateCoefficients(gasPhase.kTurb, gasPhase.epsTurb,
				*(gasPhase.turbulenceSources->turbPar));

	/*Calculation of fields for SLE matrix calculation*/
	nonIdealCorrectionOld.ref_r() = gasPhase.phaseThermodynamics->nonIdeality(
			gasPhase.concentration.p, gasPhase.temperature.ref());

	CVOld.ref_r() = gasPhase.phaseThermodynamics->Cv(gasPhase.concentration.p);

	CvM.ref_r() = CVOld.ref() * gasPhase.concentration.v[0].ref()
			/ gasPhase.density[0].ref();

	const auto surfaceCv = linearInterpolate(CVOld, bncCalc);

	CCVOld.ref_r() = gasPhase.concentration.v[0].ref() * CVOld.ref();

	const auto surfaceTemperature = linearInterpolate(gasPhase.temperature,
			bncCalc);

	concentrationsPack<quadraticSurface> surfaceConcentration { mesh_,
			gasPhase.phaseThermodynamics->Mv().size() };

	for (std::size_t k = 1; k < surfaceConcentration.v.size(); ++k)
	{
		surfaceConcentration.v[k] = linearInterpolate(
				gasPhase.concentration.v[k], bncCalc);

		surfaceConcentration.v[0].ref_r() += surfaceConcentration.v[k].ref();

		surfaceRho.ref_r() += surfaceConcentration.v[k].ref()
				* gasPhase.phaseThermodynamics->Mv()[k - 1];
	}

	const auto surfaceMav = division(surfaceRho, surfaceConcentration.v[0]);

	const auto gradMav_cell = grad(surfaceMav);

	volumeField<scalar> avMolMass { mesh_, 0 };
	{
		std::vector<boundaryConditionType> bndCon(commonConditions);
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedAverageMolarMass);
		avMolMass = volumeField<scalar>(mesh_, 0, bndCon[0], 0, bndCon[1], 0,
				bndCon[2], 0, bndCon[3], 0, bndCon[4], 0, bndCon[5], 0);
	}
	avMolMass.ref_r() = gasPhase.density[0].ref()
			/ gasPhase.concentration.v[0].ref();

	const auto gradMav = gradientLinearInterpolate(gradMav_cell, avMolMass,
			bncCalc);

	NonIdRho.ref_r() = nonIdealCorrectionOld.ref() / gasPhase.density[0].ref();

	const auto gradNonIdRho = gradientLinearInterpolate(grad(NonIdRho, bncCalc),
			NonIdRho, bncCalc);

	const auto gradCvM = gradientLinearInterpolate(grad(CvM, bncCalc), CvM,
			bncCalc);

	volumeField<tensor> gradV { mesh_ };
	volumeField<vector> gradRho { mesh_ };
	volumeField<vector> gradP { mesh_ };
	volumeField<scalar> divV { mesh_ };
	surfaceField<scalar> divV_s { mesh_ };
	volumeField<scalar> thetaS_R { mesh_ };

	if (linearRec)
	{
		gradV = grad(gasPhase.velocity, bncCalc);
		gradRho = grad(gasPhase.density[0], bncCalc);
		gradP = grad(gasPhase.pressure, bncCalc);
		divV = divergence(gasPhase.velocity, bncCalc);
	}
	else
	{
		gradV = grad(star.v);
		gradRho = grad(star.rho);
		gradP = grad(star.p);
		divV = divergence(star.v);
	}

	divV_s = divergenceLinearInterpolate(divV, diffFieldsOld.velocity, bncCalc);

	/*Calculation of all diffusion coefficients.*/
	effectiveCoeffs.pNu = harmonicInterpolate(gasPhase.pNu, bncCalc);
	effectiveCoeffs.pD = harmonicInterpolate(gasPhase.pD, bncCalc);
	effectiveCoeffs.pKappa = harmonicInterpolate(gasPhase.pKappa, bncCalc);

	if (gasPhase.turbulenceSources->turbulence)
	{
		effectiveCoeffs.tNu = harmonicInterpolate(gasPhase.tNu, bncCalc);
		effectiveCoeffs.recalculateCoefficients(
				*(gasPhase.turbulenceSources->turbPar));

		for (std::size_t k = 0; k < effectiveCoeffs.mD.size(); ++k)
			effectiveCoeffs.mD[k].ref_r() =
					gasPhase.phaseThermodynamics->Mv()[k]
							* surfaceConcentration.v[0].ref()
							* (effectiveCoeffs.pD.ref()
									+ effectiveCoeffs.tD.ref());

		effectiveCoeffs.rhoD.ref_r() = surfaceRho.ref()
				* (effectiveCoeffs.pD.ref() + effectiveCoeffs.tD.ref());

		effectiveCoeffs.mu.ref_r() = surfaceRho.ref()
				* (effectiveCoeffs.pNu.ref() + effectiveCoeffs.tNu.ref());

		effectiveCoeffs.kappa.ref_r() = effectiveCoeffs.pKappa.ref()
				+ surfaceConcentration.v[0].ref() * surfaceCv.ref()
						* effectiveCoeffs.tKappa.ref()
				+ surfaceRho.ref() * effectiveCoeffs.tLambda.ref();

		const auto k_surf = linearInterpolate(gasPhase.kTurb, bncCalc);
		const auto eps_surf = linearInterpolate(gasPhase.epsTurb, bncCalc);

		const auto thetaS_s = gasPhase.turbulenceSources->turbPar->thetaS_D(
				divV_s, k_surf, eps_surf);

		effectiveCoeffs.rhoDk.ref_r() = surfaceRho.ref()
				* (effectiveCoeffs.pNu.ref()
						+ effectiveCoeffs.k_D.ref() * thetaS_s.ref());
		effectiveCoeffs.rhoDeps.ref_r() = surfaceRho.ref()
				* (effectiveCoeffs.pNu.ref()
						+ effectiveCoeffs.eps_D.ref() * thetaS_s.ref());
		if ((gasPhase.turbulenceSources->model == turbulenceModelEnum::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::kEpsASource))
		{
			effectiveCoeffs.rhoDa.ref_r() = surfaceRho.ref() * 2.
					* (effectiveCoeffs.pNu.ref() + effectiveCoeffs.a_D.ref());
			if ((gasPhase.turbulenceSources->model
					== turbulenceModelEnum::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModelEnum::BHRKLSource))
				effectiveCoeffs.rhoDb.ref_r() =
						surfaceRho.ref()
								* (effectiveCoeffs.pNu.ref()
										+ effectiveCoeffs.b_D.ref());
		}

		thetaS_R = gasPhase.turbulenceSources->turbPar->thetaS_R(divV,
				diffFieldsOld.k, diffFieldsOld.eps);
	}
	else
	{
		for (std::size_t k = 0; k < effectiveCoeffs.mD.size(); ++k)
			effectiveCoeffs.mD[k].ref_r() =
					gasPhase.phaseThermodynamics->Mv()[k]
							* surfaceConcentration.v[0].ref()
							* effectiveCoeffs.pD.ref();

		effectiveCoeffs.rhoD.ref_r() = surfaceRho.ref()
				* effectiveCoeffs.pD.ref();

		effectiveCoeffs.mu.ref_r() = surfaceRho.ref()
				* effectiveCoeffs.pNu.ref();

		effectiveCoeffs.kappa.ref_r() = effectiveCoeffs.pKappa.ref();
	}

	if (gasPhase.turbulenceSources->turbulence)
	{
		oldField.ref_r() = diffFieldsOld.k.ref() * gasPhase.density[0].ref();
		if (msolver.solverType != matrixSolverEnum::explicitSolver)
			kMatrix.generateLaplacian(diffFieldsOld.k, oldField,
					gasPhase.density[0], effectiveCoeffs.rhoDk, timestep,
					bncCalc);
		else
			kMatrix.generateExplicitLaplacian(diffFieldsOld.k, oldField,
					gasPhase.density[0], effectiveCoeffs.rhoDk,
					linearInterpolate(gasPhase.kTurb, bncCalc), timestep,
					bncCalc);

		oldField.ref_r() = diffFieldsOld.eps.ref() * gasPhase.density[0].ref();
		if (msolver.solverType != matrixSolverEnum::explicitSolver)
			epsMatrix.generateLaplacian(diffFieldsOld.eps, oldField,
					gasPhase.density[0], effectiveCoeffs.rhoDeps, timestep,
					bncCalc);
		else
			epsMatrix.generateExplicitLaplacian(diffFieldsOld.eps, oldField,
					gasPhase.density[0], effectiveCoeffs.rhoDeps,
					linearInterpolate(gasPhase.epsTurb, bncCalc), timestep,
					bncCalc);

		if ((gasPhase.turbulenceSources->model == turbulenceModelEnum::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::kEpsASource))
		{
			if (msolver.solverType != matrixSolverEnum::explicitSolver)
				aMatrix.generateLaplacian(diffFieldsOld.a, gasPhase.density[0],
						effectiveCoeffs.rhoDa, timestep, bncCalc);
			else
			{
				if (linearRec)
					aMatrix.generateExplicitLaplacian(diffFieldsOld.a,
							gasPhase.density[0], effectiveCoeffs.rhoDa,
							linearInterpolate(diffFieldsOld.a, bncCalc),
							timestep, bncCalc);
				else
					aMatrix.generateExplicitLaplacian(diffFieldsOld.a,
							gasPhase.density[0], effectiveCoeffs.rhoDa, star.a,
							timestep, bncCalc);
			}

			if ((gasPhase.turbulenceSources->model
					== turbulenceModelEnum::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModelEnum::BHRKLSource))
			{
				oldField.ref_r() = diffFieldsOld.b.ref()
						* gasPhase.density[0].ref();
				if (msolver.solverType != matrixSolverEnum::explicitSolver)
					bMatrix.generateLaplacian(diffFieldsOld.b, oldField,
							gasPhase.density[0], effectiveCoeffs.rhoDb,
							timestep, bncCalc);
				else
				{
					if (linearRec)
						bMatrix.generateExplicitLaplacian(diffFieldsOld.b,
								oldField, gasPhase.density[0],
								effectiveCoeffs.rhoDb,
								linearInterpolate(diffFieldsOld.b, bncCalc),
								timestep, bncCalc);
					else
						bMatrix.generateExplicitLaplacian(diffFieldsOld.b,
								oldField, gasPhase.density[0],
								effectiveCoeffs.rhoDb, star.b, timestep,
								bncCalc);
				}
			}
		}
	}

	/*Generate SLE matrix for each component and calculate explicit diffusion fields.*/
	for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
	{
		oldField.ref_r() = gasPhase.density[0].ref()
				* diffFieldsOld.massFraction[k].ref();
		if (msolver.solverType != matrixSolverEnum::explicitSolver)
			massFractionMatrix[k].generateLaplacian(
					diffFieldsOld.massFraction[k], oldField,
					gasPhase.density[0], effectiveCoeffs.rhoD, timestep,
					bncCalc, k + 1);
		else
		{
			surfaceField<scalar> surfaceMassFrac_k { mesh_, 0 };

			if (linearRec)
				surfaceMassFrac_k = linearInterpolate(
						diffFieldsOld.massFraction[k], bncCalc, k + 1);
			else
				surfaceMassFrac_k.ref_r() = star.c[k + 1].ref()
						* gasPhase.phaseThermodynamics->Mv()[k]
						/ star.rho.ref();

			massFractionMatrix[k].generateExplicitLaplacian(
					diffFieldsOld.massFraction[k], oldField,
					gasPhase.density[0], effectiveCoeffs.rhoD,
					surfaceMassFrac_k, timestep, bncCalc, k + 1);
		}
	}

	std::vector<volumeField<vector>> molFracGrad(massFractionMatrix.size(),
			volumeField<vector>(mesh_, vector(0)));
	for (std::size_t k = 1; k < surfaceConcentration.v.size(); ++k)
	{
		surfaceField<scalar> surfMolFrac_k { mesh_, 0 };

		surfMolFrac_k.ref_r() = surfaceConcentration.v[k].ref()
				/ surfaceConcentration.v[0].ref();
		molFracGrad[k - 1] = grad(surfMolFrac_k);
	}

	std::vector<volumeField<scalar>> cellMolFrac(massFractionMatrix.size(),
			volumeField<scalar>(mesh_, 0));
	for (std::size_t k = 1; k < gasPhase.concentration.v.size(); ++k)
	{
		std::vector<boundaryConditionType> bndCon(commonConditions);
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedMolarFraction);
		cellMolFrac[k - 1] = volumeField<scalar>(mesh_, 0, bndCon[0], 0,
				bndCon[1], 0, bndCon[2], 0, bndCon[3], 0, bndCon[4], 0,
				bndCon[5], 0);

		cellMolFrac[k - 1].ref_r() = gasPhase.concentration.v[k].ref()
				/ gasPhase.concentration.v[0].ref();
	}

	std::vector<surfaceField<vector>> explGradW { massFractionMatrix.size(),
			surfaceField<vector>(mesh_, vector(0)) };
	std::vector<surfaceField<vector>> explGradX { massFractionMatrix.size(),
			surfaceField<vector>(mesh_, vector(0)) };

	std::vector<surfaceField<vector>> explDiffFlowX {
			gasPhase.phaseThermodynamics->Mv().size(), surfaceField<vector>(
					mesh_, vector(0)) };

	for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
	{
		explGradW[k] = gradientLinearInterpolate(
				grad(diffFieldsOld.massFraction[k], bncCalc, k + 1),
				diffFieldsOld.massFraction[k], bncCalc, k + 1);

		explGradX[k] = gradientLinearInterpolate(molFracGrad[k], cellMolFrac[k],
				bncCalc, k + 1);
	}

	for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
	{
		for (std::size_t s = 0; s < massFractionMatrix.size(); ++s)
			if (s != k)
				explDiffFlowX[k].ref_r() += astProduct(explGradX[s],
						gasPhase.phaseThermodynamics->Mv()[s]).ref();

		divisionSelf(explDiffFlowX[k], surfaceMav);
		astProductSelf(explDiffFlowX[k], effectiveCoeffs.mD[k]);

		auto explDiffFlowW_k = astProduct(effectiveCoeffs.rhoD, explGradW[k]);
		astProductSelf(explDiffFlowW_k, -1);

		/*Calculate enthalpy per mole per Kelvin (isobaric molar heat capacity) for k component*/
		surfaceField<scalar> h_k { mesh_, 0 };
		h_k.ref_r() = gasPhase.phaseThermodynamics->hkT(
				surfaceConcentration.v[k + 1].ref(), surfaceTemperature.ref(),
				k);

		/*Calculating molar mass correction for diffusion flow and explicit enthalpy flow.*/
		surfaceField<vector> massFrCorr_k { mesh_, vector(0) };
		if (molMassDiffusionFlag)
			massFrCorr_k.ref_r() = explDiffFlowX[k].ref()
					- explDiffFlowW_k.ref();

		for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
		{
			const vector entExplDiffFlow_k_i = (explDiffFlowW_k.ref()[i]
					+ massFrCorr_k.ref()[i]) * h_k.ref()[i]
					/ gasPhase.phaseThermodynamics->Mv()[k];

			entExplDiffFlow.ref_r()[i] += entExplDiffFlow_k_i;
		}

		const auto divMolMassFlow = divergence(massFrCorr_k);

		//massFractionMatrix[k].SLE[0].freeTerm -= divMolMassFlow.ref();

		massFractionMatrix[k].distributeMinusSourceTerm(divMolMassFlow,
				diffFieldsOld.massFraction[k]);
	}

	/*Calculation of resulting nonideality correction in laplacian*/
	surfaceField<vector> nonIdFlow { mesh_, vector(0) };

	nonIdFlow.ref_r() = astProduct(gradCvM, surfaceTemperature).ref()
			+ gradNonIdRho.ref();
	const auto nonIdealityCorrectionLaplacian = astProduct(nonIdFlow,
			astProduct(surfaceRho, effectiveCoeffs.tLambda));

	/*Calculation of new mass fraction*/
	{
		volumeField<scalar> sumMassFrac(mesh_, 0.);

		for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
		{
			if (msolver.solverType != matrixSolverEnum::explicitSolver)
				diffFieldsNew.massFraction[k].ref_r() = msolver.solve(
						diffFieldsOld.massFraction[k].ref(),
						massFractionMatrix[k]);
			else
				diffFieldsNew.massFraction[k].ref_r() =
						massFractionMatrix[k].SLE[0].freeTerm * timestep
								/ massFractionMatrix[k].SLE[0].centralDiagonale;

			std::replace_if(std::begin(diffFieldsNew.massFraction[k].ref_r()),
					std::end(diffFieldsNew.massFraction[k].ref_r()),
					[](const scalar value) 
					{
						return value < 0.;
					}, 0.0);

			sumMassFrac.ref_r() += diffFieldsNew.massFraction[k].ref();
		}

		for (auto & newMassFrac_k : diffFieldsNew.massFraction)
			newMassFrac_k.ref_r() /= sumMassFrac.ref();
	}

	/*Calculation of new concentrations*/
	concentrationsPack<cubicCell> intermediateConcentration { mesh_,
			gasPhase.phaseThermodynamics->Mv().size() };
	for (std::size_t k = 1; k < intermediateConcentration.v.size(); ++k)
	{
		intermediateConcentration.v[k].ref_r() = diffFieldsNew.massFraction[k
				- 1].ref() * gasPhase.density[0].ref()
				/ gasPhase.phaseThermodynamics->Mv()[k - 1];

		intermediateConcentration.v[0].ref_r() +=
				intermediateConcentration.v[k].ref();
	}

	CCVNew.ref_r() = gasPhase.phaseThermodynamics->Cv(
			intermediateConcentration.p) * intermediateConcentration.v[0].ref();
	nonIdealCorrectionNew.ref_r() = gasPhase.phaseThermodynamics->nonIdeality(
			intermediateConcentration.p, gasPhase.temperature.ref());

	oldField.ref_r() = CCVOld.ref() * diffFieldsOld.temperature.ref()
			+ nonIdealCorrectionOld.ref() - nonIdealCorrectionNew.ref();

	/*Generate SLE matrix for components of the velocity vector.*/
	surfaceField<tensor> devPhysViscSurf { mesh_, tensor(0) };

	if (msolver.solverType != matrixSolverEnum::explicitSolver)
		velocityMatrix.generateLaplacian(diffFieldsOld.velocity,
				gasPhase.density[0], effectiveCoeffs.mu, timestep, bncCalc);
	else
	{
		if (linearRec)
			velocityMatrix.generateExplicitLaplacian(diffFieldsOld.velocity,
					gasPhase.density[0], effectiveCoeffs.mu,
					linearInterpolate(diffFieldsOld.velocity, bncCalc),
					timestep, bncCalc);
		else
			velocityMatrix.generateExplicitLaplacian(diffFieldsOld.velocity,
					gasPhase.density[0], effectiveCoeffs.mu, star.v, timestep,
					bncCalc);
	}

	if (msolver.solverType != matrixSolverEnum::explicitSolver)
		temperatureMatrix.generateLaplacian(diffFieldsOld.temperature, oldField,
				CCVNew, effectiveCoeffs.kappa, timestep, bncCalc);
	else
		temperatureMatrix.generateExplicitLaplacian(diffFieldsOld.temperature,
				oldField, CCVNew, effectiveCoeffs.kappa, surfaceTemperature,
				timestep, bncCalc);

	/*Deviatoric part of physical viscosity tensor, begins from velocity gradient.*/
	const auto gradV_s = gradientLinearInterpolate(gradV, gasPhase.velocity,
			bncCalc);
	surfaceField<tensor> devTotVisc1 = gradV_s;

	volumeField<tensor> grada { mesh_, tensor(0) };
	volumeField<scalar> diva { mesh_, 0 };
	volumeField<vector> gradb { mesh_, vector(0) };

	//volumeField<vector> divDeva1 { mesh_, vector(0) };

	if (gasPhase.turbulenceSources->turbulence)
	{
		switch (gasPhase.turbulenceSources->model)
		{
		case turbulenceModelEnum::kEpsASource:
		{
			volumeField<scalar> bCalculated { mesh_, 0 };
			volumeField<scalar> numerator { mesh_, 0 };
			volumeField<scalar> denomenator { mesh_, 0 };
			for (std::size_t k = 1; k < gasPhase.concentration.v.size(); ++k)
			{
				std::valarray<scalar> x_k = gasPhase.concentration.v[k].ref()
						/ gasPhase.concentration.v[0].ref();

				numerator.ref_r() += x_k
						/ (gasPhase.density[k].ref() + stabilizator
								+ gasPhase.turbulenceSources->turbPar->Cb1()
										* gasPhase.density[0].ref());

				denomenator.ref_r() += x_k * gasPhase.density[k].ref()
						/ (gasPhase.density[k].ref() + stabilizator
								+ gasPhase.turbulenceSources->turbPar->Cb1()
										* gasPhase.density[0].ref());
			}

			bCalculated.ref_r() = gasPhase.density[0].ref() * numerator.ref()
					/ denomenator.ref() - 1.;

			std::replace_if(std::begin(bCalculated.ref_r()),
					std::end(bCalculated.ref_r()), [](const scalar value) 
					{
						return value < 0.;
					}, 0.0);

			diffFieldsOld.b.ref_r() = bCalculated.ref();
			diffFieldsNew.b.ref_r() = bCalculated.ref();

			if (linearRec)
			{
				grada = grad(gasPhase.aTurb, bncCalc);
				diva = divergence(gasPhase.aTurb, bncCalc);
			}
			else
			{
				grada = grad(star.a);
				diva = divergence(star.a);
			}

			//auto aDiffSurface = gradientLinearInterpolate(grada,
			//		diffFieldsOld.a, bncCalc);

			//for (auto & a_i : aDiffSurface.ref_r())
			//	a_i.transpose();

			//aDiffSurface.ref_r() = astProduct(aDiffSurface,
			//		effectiveCoeffs.rhoDa).ref();

			//divDeva1 = divergence(aDiffSurface);
		}
			break;
		case turbulenceModelEnum::BHRSource:
		case turbulenceModelEnum::BHRKLSource:
		{
			if (linearRec)
			{
				grada = grad(gasPhase.aTurb, bncCalc);
				diva = divergence(gasPhase.aTurb, bncCalc);
				gradb = grad(diffFieldsOld.b, bncCalc);
			}
			else
			{
				grada = grad(star.a);
				diva = divergence(star.a);
				gradb = grad(star.b);
			}

			//auto aDiffSurface = gradientLinearInterpolate(grada,
			//		diffFieldsOld.a, bncCalc);

			//for (auto & a_i : aDiffSurface.ref_r())
			//	a_i.transpose();

			//aDiffSurface.ref_r() = astProduct(aDiffSurface,
			//		effectiveCoeffs.rhoDa).ref();

			//divDeva1 = divergence(aDiffSurface);
		}
			break;
		case turbulenceModelEnum::arithmeticA1Source:
		{
			volumeField<scalar> sonicSpeed2 { mesh_, 0 };
			sonicSpeed2.ref_r() = gasPhase.phaseThermodynamics->sqSonicSpeed(
					gasPhase.concentration.p, gasPhase.density[0].ref(),
					gasPhase.internalEnergy.ref(), gasPhase.pressure.ref());

			for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
				diffFieldsOld.a.ref_r()[i] = -gasPhase.tNu.ref()[i]
						/ gasPhase.turbulenceSources->turbPar->Ca1()
						* (gradRho.ref()[i] / gasPhase.density[0].ref()[i]
								- gradP.ref()[i]
										/ (gasPhase.density[0].ref()[i]
												* sonicSpeed2.ref()[i]));

			grada = grad(diffFieldsOld.a, bncCalc);
			diva = divergence(diffFieldsOld.a, bncCalc);
		}
			break;
		case turbulenceModelEnum::arithmeticA2Source:
		{
			for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
				diffFieldsOld.a.ref_r()[i] = -gasPhase.tNu.ref()[i]
						/ gasPhase.turbulenceSources->turbPar->Ca1()
						* (gradRho.ref()[i] / gasPhase.density[0].ref()[i]
								- gradP.ref()[i] / gasPhase.pressure.ref()[i]);

			grada = grad(diffFieldsOld.a, bncCalc);
			diva = divergence(diffFieldsOld.a, bncCalc);
		}
			break;
		case turbulenceModelEnum::arithmeticA3Source:
		{
			for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
				diffFieldsOld.a.ref_r()[i] = -gasPhase.tNu.ref()[i]
						/ gasPhase.turbulenceSources->turbPar->Ca1()
						* gradRho.ref()[i] / gasPhase.density[0].ref()[i];

			grada = grad(diffFieldsOld.a, bncCalc);
			diva = divergence(diffFieldsOld.a, bncCalc);
		}
			break;
		default:
			break;
		}

		devPhysViscSurf = devTotVisc1;
	}

	/*Transposition of velocity gradient.*/
	for (auto & devTotVisc1_i : devTotVisc1.ref_r())
		devTotVisc1_i.transpose();

	if (gasPhase.turbulenceSources->turbulence)
		devPhysViscSurf.ref_r() += devTotVisc1.ref();

	/*Explicit deviatoric part of molecular viscosity tensor for energy computation.*/
	auto devPhysVisc = gradV;

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		devPhysVisc.ref_r()[i].transpose();

	devPhysVisc.ref_r() += gradV.ref();

	/*Add buoyant viscosity component to diagonal of all three tensors.*/
	for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
	{
		const scalar lambdadivV_s_i { twothirds * divV_s.ref()[i] };

		devTotVisc1.ref_r()[i].v_r()[0] -= lambdadivV_s_i;
		devTotVisc1.ref_r()[i].v_r()[4] -= lambdadivV_s_i;
		devTotVisc1.ref_r()[i].v_r()[8] -= lambdadivV_s_i;

		if (gasPhase.turbulenceSources->turbulence)
		{
			devPhysViscSurf.ref_r()[i].v_r()[0] -= lambdadivV_s_i;
			devPhysViscSurf.ref_r()[i].v_r()[4] -= lambdadivV_s_i;
			devPhysViscSurf.ref_r()[i].v_r()[8] -= lambdadivV_s_i;
		}
	}

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar lambdadivV_i { twothirds * divV.ref()[i] };

		devPhysVisc.ref_r()[i].v_r()[0] -= lambdadivV_i;
		devPhysVisc.ref_r()[i].v_r()[4] -= lambdadivV_i;
		devPhysVisc.ref_r()[i].v_r()[8] -= lambdadivV_i;
	}

	volumeField<tensor> spherTurbR { mesh_, tensor(0) };
	volumeField<tensor> devTurbR { mesh_, tensor(0) };
	volumeField<tensor> philtTurbVisc { mesh_, tensor(0) };

	if (gasPhase.turbulenceSources->turbulence)
	{
		devTurbR = devPhysVisc;
		philtTurbVisc = devPhysVisc;
	}

	/*Multiplying on physical, effective or turbulent viscosity correspondingly.*/
	devPhysVisc.ref_r() = astProduct(
			astProduct(devPhysVisc, gasPhase.density[0]), gasPhase.pNu).ref();

	if (gasPhase.turbulenceSources->turbulence)
	{
		auto turbViscCoeff1 = gasPhase.tNu;
		turbViscCoeff1.ref_r() = (1. - thetaS_R.ref()) * turbViscCoeff1.ref();

		philtTurbVisc = astProduct(
				astProduct(turbViscCoeff1, gasPhase.density[0]), philtTurbVisc);

		devTurbR = astProduct(astProduct(devTurbR, gasPhase.density[0]),
				astProduct(gasPhase.tNu, thetaS_R));

		for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		{
			/*Calculate turbulent pressure part*/
			const scalar turbPressure { -twothirds * gasPhase.rhokTurb.ref()[i] };
			spherTurbR.ref_r()[i] = tensor(turbPressure, 0, 0, 0, turbPressure,
					0, 0, 0, turbPressure);
		}
	}

	astProductSelf(devTotVisc1, effectiveCoeffs.mu);

	astProductSelf(devPhysViscSurf,
			astProduct(surfaceRho, effectiveCoeffs.pNu));

	const auto divDevPhysVisc = divergence(devPhysViscSurf);

	/*Accounting of explicit part of effective viscosity tensor*/
	const auto divDevTotVisc1 = divergence(devTotVisc1);
	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		for (std::size_t j = 0; j < vector::vsize; ++j)
			velocityMatrix.SLE[j].freeTerm[i] += divDevTotVisc1.ref()[i].v()[j];

	const auto divNonIdealityCorrectionLaplacian = divergence(
			nonIdealityCorrectionLaplacian);

	//temperatureMatrix.SLE[0].freeTerm -=
	//		divNonIdealityCorrectionLaplacian.ref();

	temperatureMatrix.distributeMinusSourceTerm(
			divNonIdealityCorrectionLaplacian, diffFieldsOld.temperature);

	if (enthalpySolverFlag == enthalpyFlowEnum::explicitSolve)
	{
		astProductSelf(entExplDiffFlow, surfaceTemperature);

		const auto divEnthFlowLaplacian = divergence(entExplDiffFlow);

		//temperatureMatrix.SLE[0].freeTerm -= divEnthFlowLaplacian.ref();

		temperatureMatrix.distributeMinusSourceTerm(divEnthFlowLaplacian,
				diffFieldsOld.temperature);
	}

	volumeField<vector> gradMav_Mav { mesh_, vector { 0 } };

	if (gasPhase.turbulenceSources->turbulence)
	{
		if (gasPhase.turbulenceSources->model
				!= turbulenceModelEnum::zeroSource)
			temperatureMatrix.SLE[0].freeTerm +=
					gasPhase.turbulenceSources->turbPar->rhoepsilon(gasPhase);

		if ((gasPhase.turbulenceSources->model == turbulenceModelEnum::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::kEpsASource))
			gradMav_Mav.ref_r() = division(gradMav_cell, avMolMass).ref();

		switch (gasPhase.turbulenceSources->model)
		{
		case turbulenceModelEnum::kEpsASource:
		case turbulenceModelEnum::BHRSource:
		case turbulenceModelEnum::BHRKLSource:
		case turbulenceModelEnum::arithmeticA1Source:
		case turbulenceModelEnum::arithmeticA2Source:
		case turbulenceModelEnum::arithmeticA3Source:
		{
			auto pDiva = astProduct(gasPhase.pressure, diva);

			pDiva.ref_r() -= dampProduct(devPhysVisc, grada).ref();

			//temperatureMatrix.SLE[0].freeTerm += pDiva.ref();

			temperatureMatrix.distributePlusSourceTerm(pDiva,
					diffFieldsOld.temperature);
		}
			break;
		default:
			break;
		}

		//FIXME Do not work correctly with MPI.
		//if ((source.model == BHRSource) || (source.model == BHRKLSource)
		//		|| (source.model == kEpsASource))
		//	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		//		for (std::size_t j = 0; j < vector::vsize; ++j)
		//			aMatrix.SLE[j].freeTerm[i] += divDeva1.ref()[i].v()[j];
	}

	temperatureMatrix.SLE[0].freeTerm += dampProduct(devPhysVisc, gradV).ref();

	if (gasPhase.turbulenceSources->turbulence)
		temperatureMatrix.SLE[0].freeTerm +=
				dampProduct(philtTurbVisc, gradV).ref();

	/*Calculating turbulent sources.*/
	volumeField<scalar> sigmaSourcek { mesh_, scalar { 0 } };
	volumeField<scalar> sigmaSourceeps { mesh_, 0 };
	volumeField<vector> sigmaSourcea { mesh_, vector(0) };
	volumeField<scalar> sigmaSourceb { mesh_, 0 };

	if (gasPhase.turbulenceSources->turbulence)
	{
		auto & nonConstMesh = const_cast<mesh&>(mesh_);

		std::tie(sigmaSourcek, sigmaSourceeps, sigmaSourcea, sigmaSourceb) =
				gasPhase.turbulenceSources->calculate(
						nonConstMesh.timestepSourceRef(), timestepCoeffs.first,
						gasPhase, diffFieldsOld, gradV, divDevPhysVisc, gradP,
						gradRho, grada, diva, gradb, spherTurbR, devTurbR,
						gradMav_Mav, *(gasPhase.phaseThermodynamics),
						gasPhase.tNu);
	}

	/*Calculation of new velocity by components of vector*/
	if (msolver.solverType != matrixSolverEnum::explicitSolver)
		diffFieldsNew.velocity.ref_r() = msolver.solve(
				diffFieldsOld.velocity.ref(), velocityMatrix);
	else
		for (std::size_t j = 0; j < vector::vsize; ++j)
			for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
				diffFieldsNew.velocity.ref_r()[i].v_r()[j] =
						velocityMatrix.SLE[j].freeTerm[i] * timestep
								/ velocityMatrix.SLE[j].centralDiagonale[i];

	/*Calculation of new temperature*/
	if (enthalpySolverFlag == enthalpyFlowEnum::implicitSolve)
	{
		if (msolver.solverType != matrixSolverEnum::explicitSolver)
			intermediateTemperature.ref_r() = msolver.solve(
					diffFieldsOld.temperature.ref(), temperatureMatrix);
		else
			intermediateTemperature.ref_r() = temperatureMatrix.SLE[0].freeTerm
					* timestep / temperatureMatrix.SLE[0].centralDiagonale;

		parallelism.correctBoundaryConditions(intermediateTemperature);

		nonIdealCorrectionOld.ref_r() =
				gasPhase.phaseThermodynamics->nonIdeality(
						gasPhase.concentration.p,
						intermediateTemperature.ref());

		oldField.ref_r() = CCVOld.ref() * intermediateTemperature.ref()
				+ nonIdealCorrectionOld.ref() - nonIdealCorrectionNew.ref();

		temperatureEnthMatrix.generateNabla(intermediateTemperature, oldField,
				CCVNew, entExplDiffFlow, timestep, bncCalc);

		diffFieldsNew.temperature.ref_r() = msolverForEnthalpy.solve(
				intermediateTemperature.ref(), temperatureEnthMatrix);
	}
	else
	{
		if (msolver.solverType != matrixSolverEnum::explicitSolver)
			diffFieldsNew.temperature.ref_r() = msolver.solve(
					diffFieldsOld.temperature.ref(), temperatureMatrix);
		else
			diffFieldsNew.temperature.ref_r() =
					temperatureMatrix.SLE[0].freeTerm * timestep
							/ temperatureMatrix.SLE[0].centralDiagonale;
	}

	if (gasPhase.turbulenceSources->turbulence)
	{
		/*Explicit source integration*/
		/*Calculation of new k*/
		if (msolver.solverType != matrixSolverEnum::explicitSolver)
			diffFieldsNew.k.ref_r() = msolver.solve(diffFieldsOld.k.ref(),
					kMatrix)
					+ astProduct(division(sigmaSourcek, gasPhase.density[0]),
							timestep).ref();
		else
			diffFieldsNew.k.ref_r() = kMatrix.SLE[0].freeTerm * timestep
					/ kMatrix.SLE[0].centralDiagonale
					+ sigmaSourcek.ref() / gasPhase.density[0].ref() * timestep;

		/*Calculation of new epsilon*/
		if (msolver.solverType != matrixSolverEnum::explicitSolver)
			diffFieldsNew.eps.ref_r() = msolver.solve(diffFieldsOld.eps.ref(),
					epsMatrix)
					+ sigmaSourceeps.ref() / gasPhase.density[0].ref()
							* timestep;
		else
			diffFieldsNew.eps.ref_r() = epsMatrix.SLE[0].freeTerm * timestep
					/ epsMatrix.SLE[0].centralDiagonale
					+ sigmaSourceeps.ref() / gasPhase.density[0].ref()
							* timestep;

		if ((gasPhase.turbulenceSources->model == turbulenceModelEnum::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::kEpsASource))
		{
			/*Calculation of new a by components of vector*/
			if (msolver.solverType != matrixSolverEnum::explicitSolver)
				diffFieldsNew.a.ref_r() = msolver.solve(diffFieldsOld.a.ref(),
						aMatrix);
			else
				for (std::size_t j = 0; j < vector::vsize; ++j)
					for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
						diffFieldsNew.a.ref_r()[i].v_r()[j] =
								aMatrix.SLE[j].freeTerm[i] * timestep
										/ aMatrix.SLE[j].centralDiagonale[i];

			diffFieldsNew.a.ref_r() +=
					astProduct(division(sigmaSourcea, gasPhase.density[0]),
							timestep).ref();

			if ((gasPhase.turbulenceSources->model
					== turbulenceModelEnum::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModelEnum::BHRKLSource))
			{
				/*Calculation of new b*/
				if (msolver.solverType != matrixSolverEnum::explicitSolver)
					diffFieldsNew.b.ref_r() = msolver.solve(
							diffFieldsOld.b.ref(), bMatrix)
							+ sigmaSourceb.ref() / gasPhase.density[0].ref()
									* timestep;
				else
					diffFieldsNew.b.ref_r() = bMatrix.SLE[0].freeTerm * timestep
							/ bMatrix.SLE[0].centralDiagonale
							+ sigmaSourceb.ref() / gasPhase.density[0].ref()
									* timestep;
			}
		}

		/*Check for non-negativity turbulent quantities.*/
		std::replace_if(std::begin(diffFieldsNew.k.ref_r()),
				std::end(diffFieldsNew.k.ref_r()),
				[&gasPhase](const scalar value) 
				{
					return value < gasPhase.turbulenceSources->turbPar->mink();
				}, gasPhase.turbulenceSources->turbPar->mink());

		std::replace_if(std::begin(diffFieldsNew.eps.ref_r()),
				std::end(diffFieldsNew.eps.ref_r()),
				[&gasPhase](
						const scalar value) 
						{
							return value < gasPhase.turbulenceSources->turbPar->mineps();
						}, gasPhase.turbulenceSources->turbPar->mineps());

		//XXX Non-physical bounding of vector a.
		if ((gasPhase.turbulenceSources->model == turbulenceModelEnum::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::kEpsASource))
			for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
				diffFieldsNew.a.ref_r()[i] = [&diffFieldsNew](
						const std::size_t i_l)
				{
					const auto module_a = diffFieldsNew.a.ref()[i_l].mag();
					const auto module_turb_V = std::sqrt(
							2 * diffFieldsNew.k.ref()[i_l]);

					if (module_a < module_turb_V)
						return diffFieldsNew.a.ref()[i_l];
					else
						return diffFieldsNew.a.ref()[i_l] / module_a
								* module_turb_V;
				}(i);

		if ((gasPhase.turbulenceSources->model == turbulenceModelEnum::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource))
			std::replace_if(std::begin(diffFieldsNew.b.ref_r()),
					std::end(diffFieldsNew.b.ref_r()),
					[&gasPhase](
							const scalar value) 
							{
								return value < gasPhase.turbulenceSources->turbPar->minb_value;
							}, gasPhase.turbulenceSources->turbPar->minb_value);
	}

	/*Recalculation of values in cells*/
	gasPhase.concentration.v[0].ref_r() = intermediateConcentration.v[0].ref();
	for (std::size_t k = 1; k < gasPhase.density.size(); ++k)
	{
		gasPhase.concentration.v[k].ref_r() =
				intermediateConcentration.v[k].ref();

		gasPhase.density[k].ref_r() = gasPhase.concentration.v[k].ref()
				* gasPhase.phaseThermodynamics->Mv()[k - 1];
	}

	gasPhase.velocity.ref_r() = diffFieldsNew.velocity.ref();

	gasPhase.momentum.ref_r() = astProduct(gasPhase.velocity,
			gasPhase.density[0]).ref();

	gasPhase.temperature.ref_r() = diffFieldsNew.temperature.ref();

	gasPhase.internalEnergy.ref_r() = gasPhase.phaseThermodynamics->UvcFromT(
			gasPhase.concentration.p, gasPhase.temperature.ref())
			* gasPhase.concentration.v[0].ref();

	gasPhase.pressure.ref_r() = gasPhase.phaseThermodynamics->pFromUv(
			gasPhase.concentration.p, gasPhase.internalEnergy.ref());

	if (gasPhase.turbulenceSources->turbulence)
	{
		gasPhase.kTurb.ref_r() = diffFieldsNew.k.ref();
		gasPhase.epsTurb.ref_r() = diffFieldsNew.eps.ref();
		if ((gasPhase.turbulenceSources->model == turbulenceModelEnum::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::kEpsASource))
		{
			gasPhase.aTurb.ref_r() = diffFieldsNew.a.ref();

			gasPhase.bTurb.ref_r() = diffFieldsNew.b.ref();
		}

		gasPhase.rhokTurb.ref_r() = astProduct(gasPhase.kTurb,
				gasPhase.density[0]).ref();
		gasPhase.rhoepsTurb.ref_r() = gasPhase.epsTurb.ref()
				* gasPhase.density[0].ref();
		if ((gasPhase.turbulenceSources->model == turbulenceModelEnum::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModelEnum::kEpsASource))
		{
			gasPhase.rhoaTurb.ref_r() = astProduct(gasPhase.aTurb,
					gasPhase.density[0]).ref();

			if ((gasPhase.turbulenceSources->model
					== turbulenceModelEnum::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModelEnum::BHRKLSource))
				gasPhase.rhobTurb.ref_r() = gasPhase.bTurb.ref()
						* gasPhase.density[0].ref();
		}
	}

	{
		const auto v2 = ampProduct(gasPhase.velocity, gasPhase.velocity);

		gasPhase.totalEnergy.ref_r() = gasPhase.internalEnergy.ref()
				+ gasPhase.density[0].ref() * v2.ref() * 0.5
				+ gasPhase.rhokTurb.ref();
	}

	gasPhase.HelmholtzEnergy.ref_r() = gasPhase.phaseThermodynamics->Fv(
			gasPhase.concentration.p, gasPhase.temperature.ref_r());

	gasPhase.entropy.ref_r() = gasPhase.phaseThermodynamics->Sv(
			gasPhase.concentration.p, gasPhase.temperature.ref());

	/*Diffusion time-step.*/
	if (sourceTimeFlag == timestepEnum::CourantAndSourceAndDiffusionTimeStep)
	{
		const auto maxVal = effectiveCoeffs.maxValue(
				gasPhase.turbulenceSources->turbulence,
				gasPhase.turbulenceSources->model);

		const auto diffuse_dt = 0.5 * timestepCoeffs.second
				/ (linearInterpolate(maxVal).ref() * minimalLengthScale.ref()
						* minimalLengthScale.ref()).min();

		auto & nonConstMesh = const_cast<mesh&>(mesh_);

		nonConstMesh.timestepSourceRef() = std::min(
				nonConstMesh.timestepSource(), diffuse_dt);
	}

	if (gasPhase.temperature.ref().min() < 0.)
		throw exception("Negative temperature after diffusion.",
				errorsEnum::negativeTemperatureError);

	const auto DiffusionEndTime { std::chrono::high_resolution_clock::now() };
	timeForDiffusion += std::chrono::duration_cast<std::chrono::milliseconds>(
			DiffusionEndTime - DiffusionStartTime).count();
}
