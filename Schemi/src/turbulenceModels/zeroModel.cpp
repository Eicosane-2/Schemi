/*
 * zeroModel.cpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#include "zeroModel.hpp"

#include <iostream>
#include <fstream>

schemi::zeroModel::zeroModel(const mesh & meshIn, const bool turb_in) :
		kEpsModels(meshIn, turb_in, false, false, turbulenceModel::zeroModel)
{
	/*Read turbulent parameters.*/
	{
		std::ifstream turbulentParametersFile { "./set/turbulentParameters.txt" };
		if (turbulentParametersFile.is_open())
			std::cout << "./set/turbulentParameters.txt is opened."
					<< std::endl;
		else
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
		minkSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		minepsilonSet(value);

		turbulentParametersFile >> skipBuffer >> value;
		CMS_RSet(value);
		turbulentParametersFile >> skipBuffer >> value;
		CMS_DSet(value);

		turbulentParametersFile.close();
	}
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
		schemi::volumeField<schemi::scalar>> schemi::zeroModel::calculate(
		scalar & sourceTimestep, const scalar,
		const bunchOfFields<cubicCell> & cellFields, const diffusiveFields&,
		const volumeField<tensor>&, const volumeField<vector>&,
		const volumeField<vector>&, const volumeField<vector>&,
		const volumeField<tensor>&, const volumeField<scalar>&,
		const volumeField<vector>&, const volumeField<tensor>&,
		const volumeField<tensor>&, const volumeField<vector>&,
		const abstractMixtureThermodynamics&,
		const volumeField<scalar>&) const noexcept
{
	auto & mesh_ { cellFields.pressure.meshRef() };

	const std::pair<volumeField<scalar>, volumeField<scalar>> Sourcek {
			volumeField<scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	const std::pair<volumeField<scalar>, volumeField<scalar>> Sourceeps {
			volumeField<scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	const std::pair<volumeField<vector>, volumeField<vector>> Sourcea {
			volumeField<vector>(mesh_, vector(0)), volumeField<vector>(mesh_,
					vector(0)) };
	const std::pair<volumeField<scalar>, volumeField<scalar>> Sourceb {
			volumeField<scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	const volumeField<scalar> gravGenField(mesh_, 0);

	sourceTimestep = std::min(mesh_.timestepSource(), veryBig);

	return std::make_tuple(Sourcek, Sourceeps, Sourcea, Sourceb, gravGenField);
}
