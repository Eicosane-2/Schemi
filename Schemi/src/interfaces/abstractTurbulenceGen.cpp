/*
 * abstractTurbulenceGen.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractTurbulenceGen.hpp"

#include <iostream>

#include "BHRGen.hpp"
#include "decayGen.hpp"
#include "shearGen.hpp"
#include "arithmeticAGen.hpp"
#include "kEpsAGen.hpp"
#include "BHRKLGen.hpp"
#include "zeroGen.hpp"

schemi::abstractTurbulenceGen::abstractTurbulenceGen(const bool turb_in,
		const turbulenceModel tm_in) noexcept :
		turbulence(turb_in), model(tm_in)
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
	if (turbulenceONString == "on")
		turbulenceFlag = true;
	else if (turbulenceONString == "off")
		turbulenceFlag = false;
	else
		throw exception("Unknown turbulence flag.",
				errors::initialisationError);

	turbulenceModel turbulenceModelFlag;
	if (sourceTypeString == "BHR")
		turbulenceModelFlag = turbulenceModel::BHRSource;
	else if (sourceTypeString == "zero")
		turbulenceModelFlag = turbulenceModel::zeroSource;
	else if (sourceTypeString == "decay")
		turbulenceModelFlag = turbulenceModel::decaySource;
	else if (sourceTypeString == "shear")
		turbulenceModelFlag = turbulenceModel::shearSource;
	else if (sourceTypeString == "arithmetic1")
		turbulenceModelFlag = turbulenceModel::arithmeticA1Source;
	else if (sourceTypeString == "arithmetic2")
		turbulenceModelFlag = turbulenceModel::arithmeticA2Source;
	else if (sourceTypeString == "kEpsA")
		turbulenceModelFlag = turbulenceModel::kEpsASource;
	else if (sourceTypeString == "arithmetic3")
		turbulenceModelFlag = turbulenceModel::arithmeticA3Source;
	else if (sourceTypeString == "BHRKL")
		turbulenceModelFlag = turbulenceModel::BHRKLSource;
	else
		throw exception("Unknown source generation flag.",
				errors::initialisationError);

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
