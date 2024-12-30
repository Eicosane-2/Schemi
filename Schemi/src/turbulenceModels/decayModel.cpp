/*
 * decayModel.cpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#include "decayModel.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>

schemi::decayModel::decayModel(const mesh & meshIn, const bool turb_in) :
		kEpsModels(meshIn, turb_in, true, true, turbulenceModel::decayModel)
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

		turbulentParametersFile.close();
	}
}

schemi::scalar schemi::decayModel::C2() const noexcept
{
	return modelParameters()[0];
}
schemi::scalar schemi::decayModel::Ca() const noexcept
{
	return modelParameters()[1];
}
schemi::scalar schemi::decayModel::Cb() const noexcept
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
		schemi::volumeField<schemi::scalar>> schemi::decayModel::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld, const volumeField<tensor>&,
		const volumeField<vector>&, const volumeField<vector>&,
		const volumeField<vector>&, const volumeField<tensor>&,
		const volumeField<scalar>&, const volumeField<vector>&,
		const volumeField<tensor>&, const volumeField<tensor>&,
		const volumeField<vector>&, const abstractMixtureThermodynamics&,
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
	const volumeField<scalar> gravGenField(mesh_, 0);

	std::valarray<scalar> modeps(diffFieldsOld.eps());
	const scalar maxeps { modeps.max() };
	std::replace_if(std::begin(modeps), std::end(modeps),
			[maxeps](const scalar value) 
			{
				return value < 1E-3 * maxeps;
			}, veryBig);

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps()[i] / diffFieldsOld.k()[i] };

		const scalar dissip(-cellFields.rhoepsTurb()[i]);

		Sourcek.first.r()[i] = 0;
		Sourcek.second.r()[i] = dissip / cellFields.kTurb()[i];

		Sourceeps.first.r()[i] = 0;
		Sourceeps.second.r()[i] = C2() * dissip / cellFields.kTurb()[i];

		/*Time-step calculation*/
		modeps[i] = std::abs(
				sourceTimestepCoeff * modeps[i]
						/ ((Sourceeps.first()[i]
								+ Sourceeps.second()[i]
										* cellFields.epsTurb()[i])
								/ cellFields.density[0]()[i] + stabilizator));

		Sourcea.first.r()[i] = 0.0;
		Sourcea.second.r()[i] = -cellFields.density[0]()[i] * ek * Ca();

		Sourceb.first.r()[i] = 0;
		Sourceb.second.r()[i] = -cellFields.density[0]()[i] * ek * Cb();
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(Sourcek, Sourceeps, Sourcea, Sourceb, gravGenField);
}
