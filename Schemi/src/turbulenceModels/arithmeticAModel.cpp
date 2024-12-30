/*
 * arithmeticAModel.cpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#include "arithmeticAModel.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>

schemi::arithmeticAModel::arithmeticAModel(const mesh & meshIn,
		const bool turb_in, const turbulenceModel model) :
		kEpsModels(meshIn, turb_in, false, false, model)
{
	modelParametersSet().resize(5);

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
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[3];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[4];

		turbulentParametersFile.close();
	}
}

schemi::scalar schemi::arithmeticAModel::C0() const noexcept
{
	return modelParameters()[0];
}
schemi::scalar schemi::arithmeticAModel::C1() const noexcept
{
	return modelParameters()[1];
}
schemi::scalar schemi::arithmeticAModel::C2() const noexcept
{
	return modelParameters()[2];
}
schemi::scalar schemi::arithmeticAModel::C3() const noexcept
{
	return modelParameters()[3];
}

schemi::scalar schemi::arithmeticAModel::sigmaRho() const noexcept
{
	return modelParameters()[4];
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
		schemi::volumeField<schemi::scalar>> schemi::arithmeticAModel::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld,
		const volumeField<tensor> & gradV,
		const volumeField<vector> & divDevPhysVisc,
		const volumeField<vector> & gradP, const volumeField<vector>&,
		const volumeField<tensor>&, const volumeField<scalar>&,
		const volumeField<vector>&, const volumeField<tensor> & spherR,
		const volumeField<tensor> & devR, const volumeField<vector>&,
		const abstractMixtureThermodynamics&,
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

		const scalar rhoSpherRGen = spherR()[i] && gradV()[i];

		const scalar rhoDevRGen = devR()[i] && gradV()[i];

		const scalar gravGen(
				diffFieldsOld.a()[i] & (gradP()[i] - divDevPhysVisc()[i]));

		const scalar dissip(-cellFields.rhoepsTurb()[i]);

		Sourcek.first.r()[i] = rhoSpherRGen + rhoDevRGen + gravGen;
		Sourcek.second.r()[i] = dissip / cellFields.kTurb()[i];

		Sourceeps.first.r()[i] = (C1() * rhoDevRGen + C3() * rhoSpherRGen
				+ C0() * gravGen) * ek;
		Sourceeps.second.r()[i] = C2() * dissip / cellFields.kTurb()[i];

		/*Time-step calculation*/
		modeps[i] = std::abs(
				sourceTimestepCoeff * modeps[i]
						/ ((Sourceeps.first()[i]
								+ Sourceeps.second()[i]
										* cellFields.epsTurb()[i])
								/ cellFields.density[0]()[i] + stabilizator));

		gravGenField.r()[i] = gravGen;
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(Sourcek, Sourceeps, Sourcea, Sourceb, gravGenField);
}

schemi::volumeField<schemi::vector> schemi::arithmeticAModel::calculate_a(
		const turbulenceModel model, const mesh & msh,
		const homogeneousPhase<cubicCell> & gasPhase,
		const volumeField<vector> & gradRho,
		const volumeField<vector> & gradP) const noexcept
{
	volumeField<vector> a(msh);

	switch (model)
	{
	case turbulenceModel::arithmeticA1Model:
	{
		volumeField<scalar> sonicSpeed2 { msh, 0 };
		sonicSpeed2.r() = gasPhase.phaseThermodynamics->sqSonicSpeed(
				gasPhase.concentration.p, gasPhase.density[0](),
				gasPhase.internalEnergy(), gasPhase.pressure());

		for (std::size_t i = 0; i < msh.cellsSize(); ++i)
			a.r()[i] = -gasPhase.tNu()[i] / sigmaRho()
					* (gradRho()[i] / gasPhase.density[0]()[i]
							- gradP()[i]
									/ (gasPhase.density[0]()[i]
											* sonicSpeed2()[i]));
	}
		break;
	case turbulenceModel::arithmeticA2Model:
	{
		for (std::size_t i = 0; i < msh.cellsSize(); ++i)
			a.r()[i] = -gasPhase.tNu()[i] / sigmaRho()
					* (gradRho()[i] / gasPhase.density[0]()[i]
							- gradP()[i] / gasPhase.pressure()[i]);
	}
		break;
	case turbulenceModel::arithmeticA3Model:
	{
		for (std::size_t i = 0; i < msh.cellsSize(); ++i)
			a.r()[i] = -gasPhase.tNu()[i] / sigmaRho() * gradRho()[i]
					/ gasPhase.density[0]()[i];
	}
		break;
	default:
		break;
	}

	return a;
}
