/*
 * BHRModel.cpp
 *
 *  Created on: 2025/03/25
 *      Author: Maxim Boldyrev
 */

#include "BHRModel.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>

#include "doubleDotProduct.hpp"
#include "intExpPow.hpp"

schemi::BHRModel::BHRModel(const mesh & meshIn, const MPIHandler & parIn,
		const volumeField<vector> & uCellIn,
		const surfaceField<vector> & uSurfIn,
		const std::pair<std::size_t, std::string> & readDataPoint,
		const bool turb_in) :
		kEpsModels(meshIn, turb_in, true, true, turbulenceModel::BHRModel), initialisation(
				meshIn, parIn, uCellIn, uSurfIn, readDataPoint)
{
	modelParametersSet().resize(12);

	/*Read turbulent parameters.*/
	{
		std::ifstream turbulentParametersFile { "./set/turbulentParameters.txt" };
		if (turbulentParametersFile.is_open())
			std::cout << "./set/turbulentParameters.txt is opened."
					<< std::endl;
		else
			[[unlikely]]
			throw std::ifstream::failure(
					"./set/turbulentParameters.txt not found.");

		std::string skipBuffer;
		scalar value;

		turbulentParametersFile >> skipBuffer >> value;
		CmuSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		sigmaScSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		sigmaTSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		sigmaESet(value);
		turbulentParametersFile >> skipBuffer >> value;
		sigmaKSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		sigmaEpsSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		sigmaaSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		sigmabSet(value);

		turbulentParametersFile >> skipBuffer >> value;
		CMS_RSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		CMS_DSet(value);

		turbulentParametersFile >> skipBuffer >> value;
		minkSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		minepsilonSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		minbSet(value);

		turbulentParametersFile >> skipBuffer >> modelParametersSet()[0];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[1];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[2];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[3];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[4];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[5];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[6];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[7];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[8];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[9];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[10];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[11];

		turbulentParametersFile.close();
	}
}

schemi::scalar schemi::BHRModel::C0() const noexcept
{
	return modelParameters()[0];
}
schemi::scalar schemi::BHRModel::C1() const noexcept
{
	return modelParameters()[1];
}
schemi::scalar schemi::BHRModel::C2() const noexcept
{
	return modelParameters()[2];
}
schemi::scalar schemi::BHRModel::C3() const noexcept
{
	return modelParameters()[3];
}
schemi::scalar schemi::BHRModel::C4() const noexcept
{
	return modelParameters()[4];
}

schemi::scalar schemi::BHRModel::Ca() const noexcept
{
	return modelParameters()[5];
}
schemi::scalar schemi::BHRModel::Cb1() const noexcept
{
	return modelParameters()[6];
}
schemi::scalar schemi::BHRModel::Cb2() const noexcept
{
	return modelParameters()[7];
}

schemi::scalar schemi::BHRModel::alpha2() const noexcept
{
	return modelParameters()[8];
}
schemi::scalar schemi::BHRModel::alpha3() const noexcept
{
	return modelParameters()[9];
}
schemi::scalar schemi::BHRModel::alpha4() const noexcept
{
	return modelParameters()[10];
}

schemi::scalar schemi::BHRModel::CMSM() const noexcept
{
	return modelParameters()[11];
}

std::tuple<
		std::pair<schemi::volumeField<schemi::scalar>,
				schemi::volumeField<schemi::scalar>>,
		std::pair<schemi::volumeField<schemi::scalar>,
				schemi::volumeField<schemi::scalar>>,
		std::pair<schemi::volumeField<schemi::vector>,
				schemi::volumeField<schemi::vector>>,
		std::pair<schemi::volumeField<schemi::scalar>,
				schemi::volumeField<schemi::scalar>>,
		schemi::volumeField<schemi::scalar>> schemi::BHRModel::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld,
		const volumeField<tensor> & gradV,
		const volumeField<vector> & divDevPhysVisc,
		const volumeField<vector> & gradP, const volumeField<vector> & gradRho,
		const volumeField<tensor> & grada, const volumeField<scalar> & diva,
		const volumeField<vector>&, const volumeField<tensor> & spherR,
		const volumeField<tensor> & devR, const volumeField<vector>&,
		const abstractMixtureThermodynamics & mixture,
		const volumeField<scalar>&) const noexcept
{
	auto & mesh_ { cellFields.pressure.meshRef() };

	std::pair<volumeField<scalar>, volumeField<scalar>> Sourcek { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	std::pair<volumeField<scalar>, volumeField<scalar>> Sourceeps { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	std::pair<volumeField<vector>, volumeField<vector>> Sourcea { volumeField<
			vector>(mesh_, vector(0)), volumeField<vector>(mesh_, vector(0)) };
	std::pair<volumeField<scalar>, volumeField<scalar>> Sourceb { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	volumeField<scalar> gravGenField(mesh_, 0);

	std::valarray<scalar> modeps(diffFieldsOld.eps.cval());
	const scalar maxeps { modeps.max() };
	std::replace_if(std::begin(modeps), std::end(modeps),
			[maxeps](const scalar value) 
			{
				return value < 1E-3 * maxeps;
			}, veryBig);

	std::valarray<scalar> mode = std::abs(cellFields.internalEnergy.cval());

	const auto a_s2 = mixture.sqSonicSpeed(cellFields.concentration.p,
			cellFields.density[0].cval(), cellFields.internalEnergy.cval(),
			cellFields.pressure.cval());

	const auto buoyancyGenerationFlag(
			initialisation.generationArea(cellFields.concentration));

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps.cval()[i]
				/ diffFieldsOld.k.cval()[i] };

		const scalar Mat2 = 2 * diffFieldsOld.k.cval()[i] / a_s2[i];
		const scalar Mat = std::sqrt(Mat2);

		const auto thetaS_i = thetaS_R(gradV.cval()[i].trace(),
				diffFieldsOld.k.cval()[i], diffFieldsOld.eps.cval()[i]);

		const scalar rhoSpherRGen(spherR.cval()[i] && gradV.cval()[i]);

		const scalar rhoDevRGen(thetaS_i * devR.cval()[i] && gradV.cval()[i]);

		const scalar gravGen(
				(diffFieldsOld.a.cval()[i]
						& (gradP.cval()[i] - divDevPhysVisc.cval()[i]))
						* buoyancyGenerationFlag.cval()[i]);

		const scalar Pi_K = (Mat * alpha3() * cellFields.rhoepsTurb.cval()[i]
				- (alpha2() * rhoDevRGen + 8 * Mat * alpha4() * rhoSpherRGen))
				* Mat;

		const scalar dissip(-cellFields.rhoepsTurb.cval()[i]);

		Sourcek.first.val()[i] = rhoSpherRGen + rhoDevRGen + gravGen + Pi_K;
		Sourcek.second.val()[i] = (1 + CMSM() * Mat2) * dissip
				/ cellFields.kTurb.cval()[i];

		Sourceeps.first.val()[i] = (C1() * rhoDevRGen + C3() * rhoSpherRGen
				+ C0() * std::max(gravGen, 0.0) + C4() * Pi_K) * ek;
		Sourceeps.second.val()[i] = C2() * dissip / cellFields.kTurb.cval()[i];

		/*Time-step calculation*/
		modeps[i] =
				std::abs(
						sourceTimestepCoeff * modeps[i]
								/ ((Sourceeps.first.cval()[i]
										+ Sourceeps.second.cval()[i]
												* cellFields.epsTurb.cval()[i])
										/ cellFields.density[0].cval()[i]
										+ stabilizator));
		mode[i] = std::abs(
				sourceTimestepCoeff * mode[i] / (gravGen + stabilizator));

		const vector bGradP(
				(gradP.cval()[i] - divDevPhysVisc.cval()[i])
						* diffFieldsOld.b.cval()[i]);

		const vector tauGradRho(
				(devR.cval()[i] * thetaS_i + spherR.cval()[i])
						/ cellFields.density[0].cval()[i] & gradRho.cval()[i]);

		const vector rhoAgradV(
				cellFields.rhoaTurb.cval()[i]
						& (grada.cval()[i] - gradV.cval()[i]));

		//const vector redistribution_a(
		//		cellFields.density[0]()[i]
		//				* ((diffFieldsOld.a()[i] & grada()[i])
		//						+ diffFieldsOld.a()[i] * diva()[i]));

		Sourcea.first.val()[i] = bGradP + tauGradRho + rhoAgradV;
		Sourcea.second.val()[i] = -cellFields.density[0].cval()[i] * ek * Ca();

		const scalar rhobDiva(-cellFields.rhobTurb.cval()[i] * diva.cval()[i]);

		const scalar baGradRho(
				-(diffFieldsOld.b.cval()[i] + 2.)
						* (diffFieldsOld.a.cval()[i] & gradRho.cval()[i]));

		//const scalar redistribution_b(cellFields.rhoaTurb()[i] & gradb()[i]);

		Sourceb.first.val()[i] = rhobDiva + baGradRho;
		Sourceb.second.val()[i] = -cellFields.density[0].cval()[i] * ek
				* (Cb1() + diffFieldsOld.b.cval()[i] * Cb2());

		gravGenField.val()[i] = gravGen;
	}

	sourceTimestep = std::min( { mesh_.timestepSource(), modeps.min(),
			mode.min() });

	return std::make_tuple(Sourcek, Sourceeps, Sourcea, Sourceb, gravGenField);
}

std::valarray<schemi::scalar> schemi::BHRModel::rhoepsilon(
		const bunchOfFields<cubicCell> & cf,
		const abstractMixtureThermodynamics & th, const volumeField<scalar> & k,
		const volumeField<scalar> & eps) const noexcept
{
	const auto a_s2 = th.sqSonicSpeed(cf.concentration.p, cf.density[0].cval(),
			cf.internalEnergy.cval(), cf.pressure.cval());

	const std::valarray<scalar> Mat2(CMSM() * 2 * k.cval() / a_s2);

	return cf.density[0].cval() * eps.cval() * (1 + Mat2);
}

void schemi::BHRModel::particlesTimeIntegration(
		const volumeField<vector> & gradRhoCell,
		const surfaceField<vector> & gradRhoSurf,
		const volumeField<vector> & uCell, const surfaceField<vector> & uSurf,
		const concentrationsPack<cubicCell> & concentrations,
		const std::vector<volumeField<scalar>> & densities,
		const boundaryConditionValue & boundVal,
		const std::valarray<scalar> & M, const scalar timestep,
		const volumeField<vector> & gradP, const volumeField<scalar> & divU,
		const volumeField<tensor> & gradU)
{
	initialisation.timeIntegration(gradRhoCell, gradRhoSurf, uCell, uSurf,
			concentrations, densities, boundVal, M, timestep, gradP, divU,
			gradU);
}

void schemi::BHRModel::particlesWriteOutput(
		const std::string & fieldDataDirectoryName, const scalar Time) const
{
	initialisation.writeOutput(fieldDataDirectoryName, Time);
}

void schemi::BHRModel::checkTransitionToTurbulenceModel(
		const volumeField<scalar> & nuCell,
		const surfaceField<scalar> & nuSurface, volumeField<scalar> & k,
		volumeField<scalar> & epsilon, volumeField<vector> & a,
		volumeField<scalar> & b,
		const concentrationsPack<cubicCell> & concentrations,
		const boundaryConditionValue & boundVal, const scalar timestep) noexcept
{
	initialisation.checkTransitionToTurbulenceModel(nuCell, nuSurface, k,
			epsilon, a, b, concentrations, boundVal, timestep);
}
