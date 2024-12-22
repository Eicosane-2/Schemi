/*
 * kEpsAModel.cpp
 *
 *  Created on: 2024/12/04
 *      Author: Maxim Boldyrev
 */

#include "kEpsAModel.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>

schemi::scalar schemi::kEpsAModel::thetaA(const vector & a, const scalar k,
		const scalar b) const noexcept
{
	if (CMSA() < aFlag)
	{
		const auto aa = (a & a);

		return aa / (2 * k) < CMSA() && b < CMSA() && aa / (2 * k) < b ?
				1 : CaMax() / Ca();
	}
	return 1;
}

schemi::kEpsAModel::kEpsAModel(const mesh & meshIn, const bool turb_in) :
		kEpsModels(meshIn, turb_in, true, false, turbulenceModel::kEpsAModel)
{
	modelParametersSet().resize(8);

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
		sigmaaSet(value);

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
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[5];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[6];
		turbulentParametersFile >> skipBuffer >> modelParametersSet()[7];

		turbulentParametersFile.close();
	}
}

schemi::scalar schemi::kEpsAModel::C0() const noexcept
{
	return modelParameters()[0];
}
schemi::scalar schemi::kEpsAModel::C1() const noexcept
{
	return modelParameters()[1];
}
schemi::scalar schemi::kEpsAModel::C2() const noexcept
{
	return modelParameters()[2];
}
schemi::scalar schemi::kEpsAModel::C3() const noexcept
{
	return modelParameters()[3];
}

schemi::scalar schemi::kEpsAModel::Ca() const noexcept
{
	return modelParameters()[4];
}

schemi::scalar schemi::kEpsAModel::CaMax() const noexcept
{
	return modelParameters()[5];
}

schemi::scalar schemi::kEpsAModel::CMSA() const noexcept
{
	return modelParameters()[6];
}

schemi::scalar schemi::kEpsAModel::chi() const noexcept
{
	return modelParameters()[7];
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
		schemi::volumeField<schemi::scalar>> schemi::kEpsAModel::calculate(
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
		const volumeField<scalar>&) const noexcept
{
	auto & mesh_ { cellFields.pressure.meshRef() };

	std::pair<volumeField<scalar>, volumeField<scalar>> Sourcek { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	std::pair<volumeField<scalar>, volumeField<scalar>> Sourceeps { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	std::pair<volumeField<vector>, volumeField<vector>> Sourcea { volumeField<
			vector>(mesh_, vector(0)), volumeField<vector>(mesh_, vector(0)) };
	const std::pair<volumeField<scalar>, volumeField<scalar>> Sourceb {
			volumeField<scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	volumeField<scalar> gravGenField(mesh_, 0);

	std::valarray<scalar> modeps(diffFieldsOld.eps());
	const scalar maxeps { modeps.max() };
	std::replace_if(std::begin(modeps), std::end(modeps),
			[maxeps](const scalar value) 
			{
				return value < 1E-3 * maxeps;
			}, veryBig);

	const auto a_s2 = mixture.sqSonicSpeed(cellFields.concentration.p,
			cellFields.density[0](), cellFields.internalEnergy(),
			cellFields.pressure());

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps()[i] / diffFieldsOld.k()[i] };

		const auto thetaS_i = thetaS_R(gradV()[i].trace(), diffFieldsOld.k()[i],
				diffFieldsOld.eps()[i]);

		const auto thetaA_i = thetaA(diffFieldsOld.a()[i], diffFieldsOld.k()[i],
				diffFieldsOld.b()[i]);

		const scalar rhoSpherRGen(spherR()[i] && gradV()[i]);

		const scalar rhoDevRGen(thetaS_i * devR()[i] && gradV()[i]);

		const scalar gravGen(
				(diffFieldsOld.a()[i] & (gradP()[i] - divDevPhysVisc()[i])));

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

		const vector bGradP(
				(gradP()[i] - divDevPhysVisc()[i]) * diffFieldsOld.b()[i]);

		const vector tauGradRho(
				(devR()[i] * thetaS_i + spherR()[i])
						/ cellFields.density[0]()[i] & gradRho()[i]);

		const vector rhoAgradV(
				cellFields.rhoaTurb()[i] & (grada()[i] - gradV()[i]));

		Sourcea.first.r()[i] = bGradP + tauGradRho + rhoAgradV;
		Sourcea.second.r()[i] = -cellFields.density[0]()[i] * ek * Ca()
				* thetaA_i;

		gravGenField.r()[i] = gravGen;
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(Sourcek, Sourceeps, Sourcea, Sourceb, gravGenField);
}

schemi::volumeField<schemi::scalar> schemi::kEpsAModel::calculate_b(
		const mesh & m, const std::vector<field<scalar, cubicCell>> & c,
		const std::vector<field<scalar, cubicCell>> & rho) const noexcept
{
	volumeField<scalar> bCalculated { m, 0 };
	volumeField<scalar> numerator { m, 0 };
	volumeField<scalar> denomenator { m, 0 };
	for (std::size_t k = 1; k < c.size(); ++k)
	{
		std::valarray<scalar> x_k = c[k]() / c[0]();

		numerator.r() += x_k / (rho[k]() + stabilizator + chi() * rho[0]());

		denomenator.r() += x_k * rho[k]()
				/ (rho[k]() + stabilizator + chi() * rho[0]());
	}

	bCalculated.r() = rho[0]() * numerator() / denomenator() - 1.;

	std::replace_if(std::begin(bCalculated.r()), std::end(bCalculated.r()),
			[](const scalar value) 
			{
				return value < 0.;
			}, 0.0);

	return bCalculated;
}
