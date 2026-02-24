/*
 * shearModel.cpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#include "shearModel.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>

#include "doubleDotProduct.hpp"

schemi::shearModel::shearModel(const mesh & meshIn, const bool turb_in) :
		kEpsModels(meshIn, turb_in, false, false, turbulenceModel::shearModel)
{
	modelParametersSet().resize(3);

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
		CMS_RSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		CMS_DSet(value);

		turbulentParametersFile >> skipBuffer >> value;
		minkSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		minepsilonSet(value);

		turbulentParametersFile >> skipBuffer >> modelParametersSet()[0];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[1];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[2];

		turbulentParametersFile.close();
	}
}

schemi::scalar schemi::shearModel::C1() const noexcept
{
	return modelParameters()[0];
}
schemi::scalar schemi::shearModel::C2() const noexcept
{
	return modelParameters()[1];
}
schemi::scalar schemi::shearModel::C3() const noexcept
{
	return modelParameters()[2];
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
		schemi::volumeField<schemi::scalar>> schemi::shearModel::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld,
		const volumeField<tensor> & gradV, const volumeField<vector>&,
		const volumeField<vector>&, const volumeField<vector>&,
		const volumeField<tensor>&, const volumeField<scalar>&,
		const volumeField<vector>&, const volumeField<tensor> & spherR,
		const volumeField<tensor> & devR, const volumeField<vector>&,
		const abstractMixtureThermodynamics&, const volumeField<scalar>&,
		const boundaryConditionValue&) const noexcept
{
	auto & mesh_ { cellFields.pressure.meshRef() };

	std::pair<volumeField<scalar>, volumeField<scalar>> Sourcek { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	std::pair<volumeField<scalar>, volumeField<scalar>> Sourceeps { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	const std::pair<volumeField<vector>, volumeField<vector>> Sourcea {
			volumeField<vector>(mesh_, vector(0)), volumeField<vector>(mesh_,
					vector(0)) };
	const std::pair<volumeField<scalar>, volumeField<scalar>> Sourceb {
			volumeField<scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	const volumeField<scalar> gravGenField(mesh_, 0);

	std::valarray<scalar> modeps(diffFieldsOld.eps.cval());
	const scalar maxeps { modeps.max() };
	std::replace_if(std::begin(modeps), std::end(modeps),
			[maxeps](const scalar value) 
			{
				return value < 1E-3 * maxeps;
			}, veryBig);

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps.cval()[i]
				/ diffFieldsOld.k.cval()[i] };

		const scalar rhoSpherRGen(spherR.cval()[i] && gradV.cval()[i]);

		const scalar rhoDevRGen(devR.cval()[i] && gradV.cval()[i]);

		const scalar dissip(-cellFields.rhoepsTurb.cval()[i]);

		Sourcek.first.val()[i] = rhoSpherRGen + rhoDevRGen;
		Sourcek.second.val()[i] = dissip / cellFields.kTurb.cval()[i];

		Sourceeps.first.val()[i] = (C1() * rhoDevRGen + C3() * rhoSpherRGen)
				* ek;
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
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(Sourcek, Sourceeps, Sourcea, Sourceb, gravGenField);
}
