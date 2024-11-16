/*
 * abstractTurbulenceGen.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractTurbulenceGen.hpp"

#include <iostream>

#include "BHRGen.hpp"
#include "BHR2Gen.hpp"
#include "decayGen.hpp"
#include "shearGen.hpp"
#include "arithmeticAGen.hpp"
#include "kEpsAGen.hpp"
#include "BHRKLGen.hpp"
#include "zeroGen.hpp"

schemi::abstractTurbulenceGen::abstractTurbulenceGen(const bool turb_in,
		const turbulenceModel tm_in) noexcept :
		turbulence(turb_in), model(tm_in), aField(
				[](const turbulenceModel modelType)
				{
					return ((modelType == turbulenceModel::BHRSource)
							|| (modelType == turbulenceModel::BHRKLSource)
							|| (modelType == turbulenceModel::kEpsASource)
							|| (modelType == turbulenceModel::BHR2Source));
				}(model)), bField(
				[](const turbulenceModel modelType)
				{
					return ((modelType == turbulenceModel::BHRSource)
							|| (modelType == turbulenceModel::BHRKLSource)
							|| (modelType == turbulenceModel::BHR2Source));
				}(model))
{
}

schemi::abstractTurbulenceGen::~abstractTurbulenceGen() noexcept
{
}

std::unique_ptr<schemi::abstractTurbulenceGen> schemi::abstractTurbulenceGen::createTurbulenceModel(
		const mesh & meshIn, const std::string & turbulenceONString,
		const std::string & sourceTypeString)
{
	bool turbulenceFlag;
	try
	{
		turbulenceFlag = onOffMap.at(turbulenceONString);
	} catch (std::out_of_range&)
	{
		throw exception("Unknown turbulence flag.",
				errors::initialisationError);
	}

	std::map<std::string, turbulenceModel> turbulenceModels;
	turbulenceModels.insert( { "BHR", turbulenceModel::BHRSource });
	turbulenceModels.insert( { "zero", turbulenceModel::zeroSource });
	turbulenceModels.insert( { "decay", turbulenceModel::decaySource });
	turbulenceModels.insert( { "shear", turbulenceModel::shearSource });
	turbulenceModels.insert( { "arithmetic1",
			turbulenceModel::arithmeticA1Source });
	turbulenceModels.insert( { "arithmetic2",
			turbulenceModel::arithmeticA2Source });
	turbulenceModels.insert( { "kEpsA", turbulenceModel::kEpsASource });
	turbulenceModels.insert( { "arithmetic3",
			turbulenceModel::arithmeticA3Source });
	turbulenceModels.insert( { "BHRKL", turbulenceModel::BHRKLSource });
	turbulenceModels.insert( { "BHR2", turbulenceModel::BHR2Source });

	turbulenceModel turbulenceModelFlag;
	try
	{
		turbulenceModelFlag = turbulenceModels.at(sourceTypeString);
	} catch (std::out_of_range&)
	{
		throw exception("Unknown source generation flag.",
				errors::initialisationError);
	}

	std::unique_ptr<abstractTurbulenceGen> turbulenceSources;
	switch (turbulenceModelFlag)
	{
	case turbulenceModel::BHRSource:
		turbulenceSources = std::make_unique<BHRGen>(meshIn, turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::decaySource:
		turbulenceSources = std::make_unique<decayGen>(meshIn, turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::shearSource:
		turbulenceSources = std::make_unique<shearGen>(meshIn, turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::arithmeticA1Source:
	case turbulenceModel::arithmeticA2Source:
	case turbulenceModel::arithmeticA3Source:
		turbulenceSources = std::make_unique<arithmeticAGen>(meshIn,
				turbulenceFlag, turbulenceModelFlag);
		break;
	case turbulenceModel::kEpsASource:
		turbulenceSources = std::make_unique<kEpsAGen>(meshIn, turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::BHRKLSource:
		turbulenceSources = std::make_unique<BHRKLGen>(meshIn, turbulenceFlag,
				turbulenceModelFlag);
		break;
	case turbulenceModel::BHR2Source:
		turbulenceSources = std::make_unique<BHR2Gen>(meshIn, turbulenceFlag,
				turbulenceModelFlag);
		break;
	default:
		turbulenceSources = std::make_unique<zeroGen>(meshIn, turbulenceFlag,
				turbulenceModelFlag);
		break;
	}

	/*Read turbulent parameters.*/
	{
		std::ifstream turbulentParametersFile { "./set/turbulentParameters.txt" };
		if (turbulentParametersFile.is_open())
			std::cout << "./set/turbulentParameters.txt is opened."
					<< std::endl;
		else
			throw std::ifstream::failure(
					"./set/turbulentParameters.txt not found.");

		turbulenceSources->turbPar->readTurbulentParameters(
				turbulentParametersFile);

		turbulentParametersFile.close();
	}

	return turbulenceSources;
}
