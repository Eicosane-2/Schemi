/*
 * BHRKLModel.cpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#include "BHRKLModel.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>

#include "doubleDotProduct.hpp"

schemi::BHRKLModel::BHRKLModel(const mesh & meshIn, const bool turb_in) :
		kLModels(meshIn, turb_in, true, true, turbulenceModel::BHRKLModel)
{
	modelParametersSet().resize(6);

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

		turbulentParametersFile.close();
	}
}

schemi::scalar schemi::BHRKLModel::C0() const noexcept
{
	return modelParameters()[0];
}
schemi::scalar schemi::BHRKLModel::C1() const noexcept
{
	return modelParameters()[1];
}
schemi::scalar schemi::BHRKLModel::C2() const noexcept
{
	return modelParameters()[2];
}
schemi::scalar schemi::BHRKLModel::C3() const noexcept
{
	return modelParameters()[3];
}

schemi::scalar schemi::BHRKLModel::Ca() const noexcept
{
	return modelParameters()[4];
}
schemi::scalar schemi::BHRKLModel::Cb() const noexcept
{
	return modelParameters()[5];
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
		schemi::volumeField<schemi::scalar>> schemi::BHRKLModel::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld,
		const volumeField<tensor> & gradV,
		const volumeField<vector> & divDevPhysVisc,
		const volumeField<vector> & gradP, const volumeField<vector> & gradRho,
		const volumeField<tensor> & grada, const volumeField<scalar>&,
		const volumeField<vector>&, const volumeField<tensor> & spherR,
		const volumeField<tensor> & devR, const volumeField<vector>&,
		const abstractMixtureThermodynamics & mixture,
		const volumeField<scalar>&,
		const boundaryConditionValue&) const noexcept
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

	const auto a_s2 = mixture.sqSonicSpeed(cellFields.concentration.p,
			cellFields.density[0].cval(), cellFields.internalEnergy.cval(),
			cellFields.pressure.cval());

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps.cval()[i]
				/ diffFieldsOld.k.cval()[i] };
		const scalar sqrtke { std::sqrt(diffFieldsOld.k.cval()[i])
				/ diffFieldsOld.eps.cval()[i] };

		const auto thetaS_i = thetaS_R(gradV.cval()[i].trace(),
				diffFieldsOld.k.cval()[i], diffFieldsOld.eps.cval()[i]);

		const scalar rhoSpherRGen(spherR.cval()[i] && gradV.cval()[i]);

		const scalar rhoDevRGen(thetaS_i * devR.cval()[i] && gradV.cval()[i]);

		const scalar gravGen(
				(diffFieldsOld.a.cval()[i]
						& (gradP.cval()[i] - divDevPhysVisc.cval()[i])));

		const scalar dissip(
				-cellFields.density[0].cval()[i] * diffFieldsOld.k.cval()[i]
						* sqrtke);

		Sourcek.first.val()[i] = rhoSpherRGen + rhoDevRGen + gravGen;
		Sourcek.second.val()[i] = dissip / cellFields.kTurb.cval()[i];

		Sourceeps.first.val()[i] = (C1() * rhoDevRGen + C3() * rhoSpherRGen
				+ C0() * gravGen) * ek;
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

		const vector bGradP(
				(gradP.cval()[i] - divDevPhysVisc.cval()[i])
						* diffFieldsOld.b.cval()[i]);

		const vector tauGradRho(
				(devR.cval()[i] * thetaS_i + spherR.cval()[i])
						/ cellFields.density[0].cval()[i] & gradRho.cval()[i]);

		const vector rhoAgradV(
				cellFields.rhoaTurb.cval()[i]
						& (grada.cval()[i] - gradV.cval()[i]));

		Sourcea.first.val()[i] = bGradP + tauGradRho + rhoAgradV;
		Sourcea.second.val()[i] = -cellFields.density[0].cval()[i] * sqrtke
				* Ca();

		const scalar bagradRho(
				-2. * (diffFieldsOld.b.cval()[i] + 1.)
						* (diffFieldsOld.a.cval()[i] & gradRho.cval()[i]));

		Sourceb.first.val()[i] = bagradRho;
		Sourceb.second.val()[i] = -cellFields.density[0].cval()[i] * sqrtke
				* Cb();

		gravGenField.val()[i] = gravGen;
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(Sourcek, Sourceeps, Sourcea, Sourceb, gravGenField);
}
