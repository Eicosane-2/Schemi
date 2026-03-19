/*
 * abstractMixtureThermodynamics.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractMixtureThermodynamics.hpp"

#include <fstream>
#include <iostream>

#include "exception.hpp"
#include "gasModelEnum.hpp"
#include "globalConstants.hpp"
#include "mixtureVanDerWaals.hpp"
#include "mixtureRedlichKwong.hpp"
#include "mixtureStiffened.hpp"
#include "mixtureKataokaVanDerWaals.hpp"
#include "mixtureIdeal.hpp"

schemi::abstractMixtureThermodynamics::abstractMixtureThermodynamics(
		const scalar Rin, const scalar hPin) noexcept :
		R(Rin), kB(R / NAvogardro), hPlanck(hPin)
{
}

schemi::abstractMixtureThermodynamics::~abstractMixtureThermodynamics() noexcept
{
}

std::unique_ptr<schemi::abstractMixtureThermodynamics> schemi::abstractMixtureThermodynamics::createThermodynamics(
		const std::string & equationOfState, const scalar R,
		const scalar hPlanck,
		std::array<std::valarray<scalar>, 4> & thermodynamicalProperties,
		const std::size_t numberOfComponents)
{
	std::map<std::string, gasModel> gasModels;
	gasModels.insert( { "vanDerWaals", gasModel::vanDerWaals });
	gasModels.insert( { "ideal", gasModel::ideal });
	gasModels.insert( { "RedlichKwong", gasModel::RedlichKwong });
	gasModels.insert( { "stiffened", gasModel::stiffened });
	gasModels.insert( { "Kataoka", gasModel::KataokaVanDerWaals });
	gasModel gasModelFlag;
	try
	{
		gasModelFlag = gasModels.at(equationOfState);
	} catch (const std::out_of_range&)
	{
		throw exception("Unknown type of equation of state.",
				errors::initialisationError);
	}

	std::string skipBuffer;

	switch (gasModelFlag)
	{
	case gasModel::vanDerWaals:
		return std::make_unique<mixtureVanDerWaals>(R, hPlanck,
				std::get<0>(thermodynamicalProperties),
				std::get<1>(thermodynamicalProperties),
				std::get<2>(thermodynamicalProperties),
				std::get<3>(thermodynamicalProperties));
		break;
	case gasModel::RedlichKwong:
		return std::make_unique<mixtureRedlichKwong>(R, hPlanck,
				std::get<0>(thermodynamicalProperties),
				std::get<1>(thermodynamicalProperties),
				std::get<2>(thermodynamicalProperties),
				std::get<3>(thermodynamicalProperties));
		break;
	case gasModel::stiffened:
	{
		const std::string stiffenedCoeffsName { "./set/stiffenedFluidData.txt" };

		std::ifstream fluidConditionsFile { stiffenedCoeffsName };

		if (fluidConditionsFile.is_open())
			std::cout << stiffenedCoeffsName << " is opened." << std::endl;
		else
			[[unlikely]]
			throw std::ifstream::failure(stiffenedCoeffsName + " not found.");

		std::valarray<scalar> p0Data(numberOfComponents), gammaData(
				numberOfComponents);

		for (std::size_t k = 0; k < numberOfComponents; ++k)
		{
			if (fluidConditionsFile.eof())
				throw std::ifstream::failure(
						stiffenedCoeffsName + ". Unexpected end of file.");

			fluidConditionsFile >> skipBuffer >> p0Data[k] >> gammaData[k];
		}

		return std::make_unique<mixtureStiffened>(R, hPlanck,
				std::get<0>(thermodynamicalProperties),
				std::get<1>(thermodynamicalProperties), p0Data, gammaData);
	}
		break;
	case gasModel::KataokaVanDerWaals:
	{
		const std::string KataokaCoeffsName {
				"./set/KataokaVanDerWaalsFluidData.txt" };

		std::ifstream fluidConditionsFile { KataokaCoeffsName };

		if (fluidConditionsFile.is_open())
			std::cout << KataokaCoeffsName << " is opened." << std::endl;
		else
			[[unlikely]]
			throw std::ifstream::failure(KataokaCoeffsName + " not found.");

		std::valarray<scalar> epsilonLJData(numberOfComponents), sigmaLJData(
				numberOfComponents);

		for (std::size_t k = 0; k < numberOfComponents; ++k)
		{
			if (fluidConditionsFile.eof())
				throw std::ifstream::failure(
						KataokaCoeffsName + ". Unexpected end of file.");

			fluidConditionsFile >> skipBuffer >> epsilonLJData[k]
					>> sigmaLJData[k];
		}

		std::string bCalcTypeStr;
		scalar bCoeff;

		fluidConditionsFile >> skipBuffer >> bCalcTypeStr >> bCoeff;

		std::pair<bool, scalar> bCalc;

		if (bCalcTypeStr == "critical")
			bCalc = { true, bCoeff };
		else if (bCalcTypeStr == "Kataoka")
			bCalc = { false, bCoeff };
		else
			[[unlikely]]
			throw exception(
					"Wrong type of Kataoka-van der Waals b coefficient calculation",
					errors::initialisationError);

		return std::make_unique<mixtureKataokaVanDerWaals>(R, hPlanck,
				std::get<0>(thermodynamicalProperties),
				std::get<1>(thermodynamicalProperties),
				std::get<2>(thermodynamicalProperties),
				std::get<3>(thermodynamicalProperties), epsilonLJData,
				sigmaLJData, bCalc);
	}
		break;
	case gasModel::ideal:
	default:
		return std::make_unique<mixtureIdeal>(R, hPlanck,
				std::get<0>(thermodynamicalProperties),
				std::get<1>(thermodynamicalProperties));
		break;
	}
}
