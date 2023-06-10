/*
 * boundaryConditionFromString.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "boundaryConditionFromString.hpp"

#include <iostream>
#include <algorithm>

#include "exception.hpp"

void schemi::boundaryConditionFromString(
		const std::string & boundaryConditionString,
		const std::string & boundaryConditionValueString,
		boundaryConditionType & boundaryCondition,
		std::vector<scalar> & boundaryConditionValue)
{
	boundaryConditionValue.resize(1);

	if (boundaryConditionString == "blank")
	{
		boundaryCondition = boundaryConditionType::blank;
		boundaryConditionValue[0] = 0;
	}
	else if (boundaryConditionString == "freeBoundary")
	{
		boundaryCondition = boundaryConditionType::freeBoundary;
		boundaryConditionValue[0] = 0;
	}
	else if (boundaryConditionString == "slip")
	{
		boundaryCondition = boundaryConditionType::slip;
		boundaryConditionValue[0] = 0;
	}
	else if (boundaryConditionString == "fixedValue")
	{
		boundaryCondition = boundaryConditionType::fixedValue;

		const std::size_t numberOfBoundaryZones = std::count(
				boundaryConditionValueString.begin(),
				boundaryConditionValueString.end(), ';') + 1;

		if (numberOfBoundaryZones == 1)
		{
			boundaryConditionValue.resize(1);

			boundaryConditionValue[0] = std::stod(boundaryConditionValueString);
		}
		else
		{
			boundaryConditionValue.resize(numberOfBoundaryZones);

			auto numBeginIter = boundaryConditionValueString.begin();
			auto numEndIter = boundaryConditionValueString.end();

			for (auto striter = boundaryConditionValueString.begin();
					striter != boundaryConditionValueString.end(); ++striter)
			{
				numEndIter = striter;

				if (*striter == ';')
					break;
			}

			boundaryConditionValue[0] = std::stod(
					std::string(numBeginIter, numEndIter));

			for (std::size_t i = 1; i < numberOfBoundaryZones; ++i)
			{
				numBeginIter = numEndIter;
				numBeginIter++;

				for (auto striter = numBeginIter;
						striter != boundaryConditionValueString.end();
						++striter)
				{
					numEndIter = striter;

					if (*striter == ';')
						break;
				}

				boundaryConditionValue[i] = std::stod(
						std::string(numBeginIter, numEndIter));
			}
		}
	}
	else if (boundaryConditionString == "innerSurface")
		throw exception("<<innerSurface>> can not be boundary surface type.",
				errorsEnum::boundaryConditionError);
	else
		throw exception("Unknown type of boundary condition",
				errorsEnum::boundaryConditionError);
}

void schemi::boundaryConditionFromString(
		const std::string & boundaryConditionString,
		const std::string & boundaryConditionValueStringX,
		const std::string & boundaryConditionValueStringY,
		const std::string & boundaryConditionValueStringZ,
		boundaryConditionType & boundaryCondition,
		std::vector<vector> & boundaryConditionValue)
{
	boundaryConditionValue.resize(1);

	if (boundaryConditionString == "blank")
	{
		boundaryCondition = boundaryConditionType::blank;
		boundaryConditionValue[0] = vector(0);
	}
	else if (boundaryConditionString == "freeBoundary")
	{
		boundaryCondition = boundaryConditionType::freeBoundary;
		boundaryConditionValue[0] = vector(0);
	}
	else if (boundaryConditionString == "slip")
	{
		boundaryCondition = boundaryConditionType::slip;
		boundaryConditionValue[0] = vector(0);
	}
	else if (boundaryConditionString == "fixedValue")
	{
		boundaryCondition = boundaryConditionType::fixedValue;

		const std::size_t numberOfBoundaryZonesX = std::count(
				boundaryConditionValueStringX.begin(),
				boundaryConditionValueStringX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesY = std::count(
				boundaryConditionValueStringY.begin(),
				boundaryConditionValueStringY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZ = std::count(
				boundaryConditionValueStringZ.begin(),
				boundaryConditionValueStringZ.end(), ';') + 1;

		if (!((numberOfBoundaryZonesX == numberOfBoundaryZonesY)
				&& (numberOfBoundaryZonesY == numberOfBoundaryZonesZ)))
			throw exception("Inequal number of vector boundary conditions.",
					errorsEnum::boundaryConditionError);

		if (numberOfBoundaryZonesX == 1)
			boundaryConditionValue[0] = vector(
					std::stod(boundaryConditionValueStringX),
					std::stod(boundaryConditionValueStringY),
					std::stod(boundaryConditionValueStringZ));
		else
		{
			boundaryConditionValue.resize(numberOfBoundaryZonesX);

			std::array<std::string, 3> boundaryConditionValueStringV {
					boundaryConditionValueStringX,
					boundaryConditionValueStringY, boundaryConditionValueStringZ };

			for (std::size_t k = 0; k < vector::vsize; ++k)
			{
				auto numBeginIter_k = boundaryConditionValueStringV[k].begin();
				auto numEndIter_k = boundaryConditionValueStringV[k].end();

				for (auto striter = boundaryConditionValueStringV[k].begin();
						striter != boundaryConditionValueStringV[k].end();
						++striter)
				{
					numEndIter_k = striter;

					if (*striter == ';')
						break;
				}

				boundaryConditionValue[0].v_r()[k] = std::stod(
						std::string(numBeginIter_k, numEndIter_k));

				for (std::size_t i = 1; i < numberOfBoundaryZonesX; ++i)
				{
					numBeginIter_k = numEndIter_k;
					numBeginIter_k++;

					for (auto striter = numBeginIter_k;
							striter != boundaryConditionValueStringV[k].end();
							++striter)
					{
						numEndIter_k = striter;

						if (*striter == ';')
							break;
					}

					boundaryConditionValue[i].v_r()[k] = std::stod(
							std::string(numBeginIter_k, numEndIter_k));
				}
			}
		}
	}
	else if (boundaryConditionString == "innerSurface")
		throw exception("<<innerSurface>> can not be boundary surface type.",
				errorsEnum::boundaryConditionError);
	else
		throw exception("Unknown type of boundary condition",
				errorsEnum::boundaryConditionError);
}

void schemi::boundaryConditionFromString(
		const std::string & boundaryConditionString,
		const std::string & boundaryConditionValueStringXX,
		const std::string & boundaryConditionValueStringXY,
		const std::string & boundaryConditionValueStringXZ,
		const std::string & boundaryConditionValueStringYX,
		const std::string & boundaryConditionValueStringYY,
		const std::string & boundaryConditionValueStringYZ,
		const std::string & boundaryConditionValueStringZX,
		const std::string & boundaryConditionValueStringZY,
		const std::string & boundaryConditionValueStringZZ,
		boundaryConditionType & boundaryCondition,
		std::vector<tensor> & boundaryConditionValue)
{
	boundaryConditionValue.resize(1);

	if (boundaryConditionString == "blank")
	{
		boundaryCondition = boundaryConditionType::blank;
		boundaryConditionValue[0] = tensor(0);
	}
	else if (boundaryConditionString == "freeBoundary")
	{
		boundaryCondition = boundaryConditionType::freeBoundary;
		boundaryConditionValue[0] = tensor(0);
	}
	else if (boundaryConditionString == "slip")
	{
		boundaryCondition = boundaryConditionType::slip;
		boundaryConditionValue[0] = tensor(0);
	}
	else if (boundaryConditionString == "fixedValue")
	{
		boundaryCondition = boundaryConditionType::fixedValue;

		const std::size_t numberOfBoundaryZonesXX = std::count(
				boundaryConditionValueStringXX.begin(),
				boundaryConditionValueStringXX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXY = std::count(
				boundaryConditionValueStringXY.begin(),
				boundaryConditionValueStringXY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXZ = std::count(
				boundaryConditionValueStringXZ.begin(),
				boundaryConditionValueStringXZ.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYX = std::count(
				boundaryConditionValueStringYX.begin(),
				boundaryConditionValueStringYX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYY = std::count(
				boundaryConditionValueStringYY.begin(),
				boundaryConditionValueStringYY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYZ = std::count(
				boundaryConditionValueStringYZ.begin(),
				boundaryConditionValueStringYZ.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZX = std::count(
				boundaryConditionValueStringZX.begin(),
				boundaryConditionValueStringZX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZY = std::count(
				boundaryConditionValueStringZY.begin(),
				boundaryConditionValueStringZY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZZ = std::count(
				boundaryConditionValueStringZZ.begin(),
				boundaryConditionValueStringZZ.end(), ';') + 1;

		if (!(

		(numberOfBoundaryZonesXX == numberOfBoundaryZonesXY)
				&& (numberOfBoundaryZonesXY == numberOfBoundaryZonesXZ)
				&& (numberOfBoundaryZonesXZ == numberOfBoundaryZonesYX)
				&& (numberOfBoundaryZonesYX == numberOfBoundaryZonesYY)
				&& (numberOfBoundaryZonesYY == numberOfBoundaryZonesYZ)
				&& (numberOfBoundaryZonesYZ == numberOfBoundaryZonesZX)
				&& (numberOfBoundaryZonesZX == numberOfBoundaryZonesZY)
				&& (numberOfBoundaryZonesZY == numberOfBoundaryZonesZZ)

		))
			throw exception("Inequal number of tensor boundary conditions.",
					errorsEnum::boundaryConditionError);

		if (numberOfBoundaryZonesXX == 1)
			boundaryConditionValue[0] = tensor(
					std::stod(boundaryConditionValueStringXX),
					std::stod(boundaryConditionValueStringXY),
					std::stod(boundaryConditionValueStringXZ),
					std::stod(boundaryConditionValueStringYX),
					std::stod(boundaryConditionValueStringYY),
					std::stod(boundaryConditionValueStringYZ),
					std::stod(boundaryConditionValueStringZX),
					std::stod(boundaryConditionValueStringZY),
					std::stod(boundaryConditionValueStringZZ));
		else
		{
			boundaryConditionValue.resize(numberOfBoundaryZonesXX);

			std::array<std::string, 9> boundaryConditionValueStringV {
					boundaryConditionValueStringXX,
					boundaryConditionValueStringXY,
					boundaryConditionValueStringXZ,
					boundaryConditionValueStringYX,
					boundaryConditionValueStringYY,
					boundaryConditionValueStringYZ,
					boundaryConditionValueStringZX,
					boundaryConditionValueStringZY,
					boundaryConditionValueStringZZ };

			for (std::size_t k = 0; k < tensor::vsize; ++k)
			{
				auto numBeginIter_k = boundaryConditionValueStringV[k].begin();
				auto numEndIter_k = boundaryConditionValueStringV[k].end();

				for (auto striter = boundaryConditionValueStringV[k].begin();
						striter != boundaryConditionValueStringV[k].end();
						++striter)
				{
					numEndIter_k = striter;

					if (*striter == ';')
						break;
				}

				boundaryConditionValue[0].v_r()[k] = std::stod(
						std::string(numBeginIter_k, numEndIter_k));

				for (std::size_t i = 1; i < numberOfBoundaryZonesXX; ++i)
				{
					numBeginIter_k = numEndIter_k;
					numBeginIter_k++;

					for (auto striter = numBeginIter_k;
							striter != boundaryConditionValueStringV[k].end();
							++striter)
					{
						numEndIter_k = striter;

						if (*striter == ';')
							break;
					}

					boundaryConditionValue[i].v_r()[k] = std::stod(
							std::string(numBeginIter_k, numEndIter_k));
				}
			}
		}
	}
	else if (boundaryConditionString == "innerSurface")
		throw exception("<<innerSurface>> can not be boundary surface type.",
				errorsEnum::boundaryConditionError);
	else
		throw exception("Unknown type of boundary condition",
				errorsEnum::boundaryConditionError);
}

void schemi::boundaryConditionFromString(
		const std::string & boundaryConditionString,

		const std::string & boundaryConditionValueStringXXX,
		const std::string & boundaryConditionValueStringXXY,
		const std::string & boundaryConditionValueStringXXZ,
		const std::string & boundaryConditionValueStringXYX,
		const std::string & boundaryConditionValueStringXYY,
		const std::string & boundaryConditionValueStringXYZ,
		const std::string & boundaryConditionValueStringXZX,
		const std::string & boundaryConditionValueStringXZY,
		const std::string & boundaryConditionValueStringXZZ,

		const std::string & boundaryConditionValueStringYXX,
		const std::string & boundaryConditionValueStringYXY,
		const std::string & boundaryConditionValueStringYXZ,
		const std::string & boundaryConditionValueStringYYX,
		const std::string & boundaryConditionValueStringYYY,
		const std::string & boundaryConditionValueStringYYZ,
		const std::string & boundaryConditionValueStringYZX,
		const std::string & boundaryConditionValueStringYZY,
		const std::string & boundaryConditionValueStringYZZ,

		const std::string & boundaryConditionValueStringZXX,
		const std::string & boundaryConditionValueStringZXY,
		const std::string & boundaryConditionValueStringZXZ,
		const std::string & boundaryConditionValueStringZYX,
		const std::string & boundaryConditionValueStringZYY,
		const std::string & boundaryConditionValueStringZYZ,
		const std::string & boundaryConditionValueStringZZX,
		const std::string & boundaryConditionValueStringZZY,
		const std::string & boundaryConditionValueStringZZZ,

		boundaryConditionType & boundaryCondition,
		std::vector<tensor3> & boundaryConditionValue)
{
	boundaryConditionValue.resize(1);

	if (boundaryConditionString == "blank")
	{
		boundaryCondition = boundaryConditionType::blank;
		boundaryConditionValue[0] = tensor3(0);
	}
	else if (boundaryConditionString == "freeBoundary")
	{
		boundaryCondition = boundaryConditionType::freeBoundary;
		boundaryConditionValue[0] = tensor3(0);
	}
	else if (boundaryConditionString == "slip")
	{
		boundaryCondition = boundaryConditionType::slip;
		boundaryConditionValue[0] = tensor3(0);
	}
	else if (boundaryConditionString == "fixedValue")
	{
		boundaryCondition = boundaryConditionType::fixedValue;

		const std::size_t numberOfBoundaryZonesXXX = std::count(
				boundaryConditionValueStringXXX.begin(),
				boundaryConditionValueStringXXX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXXY = std::count(
				boundaryConditionValueStringXXY.begin(),
				boundaryConditionValueStringXXY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXXZ = std::count(
				boundaryConditionValueStringXXZ.begin(),
				boundaryConditionValueStringXXZ.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXYX = std::count(
				boundaryConditionValueStringXYX.begin(),
				boundaryConditionValueStringXYX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXYY = std::count(
				boundaryConditionValueStringXYY.begin(),
				boundaryConditionValueStringXYY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXYZ = std::count(
				boundaryConditionValueStringXYZ.begin(),
				boundaryConditionValueStringXYZ.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXZX = std::count(
				boundaryConditionValueStringXZX.begin(),
				boundaryConditionValueStringXZX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXZY = std::count(
				boundaryConditionValueStringXZY.begin(),
				boundaryConditionValueStringXZY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesXZZ = std::count(
				boundaryConditionValueStringXZZ.begin(),
				boundaryConditionValueStringXZZ.end(), ';') + 1;

		const std::size_t numberOfBoundaryZonesYXX = std::count(
				boundaryConditionValueStringYXX.begin(),
				boundaryConditionValueStringYXX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYXY = std::count(
				boundaryConditionValueStringYXY.begin(),
				boundaryConditionValueStringYXY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYXZ = std::count(
				boundaryConditionValueStringYXZ.begin(),
				boundaryConditionValueStringYXZ.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYYX = std::count(
				boundaryConditionValueStringYYX.begin(),
				boundaryConditionValueStringYYX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYYY = std::count(
				boundaryConditionValueStringYYY.begin(),
				boundaryConditionValueStringYYY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYYZ = std::count(
				boundaryConditionValueStringYYZ.begin(),
				boundaryConditionValueStringYYZ.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYZX = std::count(
				boundaryConditionValueStringYZX.begin(),
				boundaryConditionValueStringYZX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYZY = std::count(
				boundaryConditionValueStringYZY.begin(),
				boundaryConditionValueStringYZY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesYZZ = std::count(
				boundaryConditionValueStringYZZ.begin(),
				boundaryConditionValueStringYZZ.end(), ';') + 1;

		const std::size_t numberOfBoundaryZonesZXX = std::count(
				boundaryConditionValueStringYXX.begin(),
				boundaryConditionValueStringYXX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZXY = std::count(
				boundaryConditionValueStringYXY.begin(),
				boundaryConditionValueStringYXY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZXZ = std::count(
				boundaryConditionValueStringYXZ.begin(),
				boundaryConditionValueStringYXZ.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZYX = std::count(
				boundaryConditionValueStringYYX.begin(),
				boundaryConditionValueStringYYX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZYY = std::count(
				boundaryConditionValueStringYYY.begin(),
				boundaryConditionValueStringYYY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZYZ = std::count(
				boundaryConditionValueStringYYZ.begin(),
				boundaryConditionValueStringYYZ.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZZX = std::count(
				boundaryConditionValueStringYZX.begin(),
				boundaryConditionValueStringYZX.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZZY = std::count(
				boundaryConditionValueStringYZY.begin(),
				boundaryConditionValueStringYZY.end(), ';') + 1;
		const std::size_t numberOfBoundaryZonesZZZ = std::count(
				boundaryConditionValueStringYZZ.begin(),
				boundaryConditionValueStringYZZ.end(), ';') + 1;

		if (!(

		(numberOfBoundaryZonesXXX == numberOfBoundaryZonesXXY)
				&& (numberOfBoundaryZonesXXY == numberOfBoundaryZonesXXZ)
				&& (numberOfBoundaryZonesXXZ == numberOfBoundaryZonesXYX)
				&& (numberOfBoundaryZonesXYX == numberOfBoundaryZonesXYY)
				&& (numberOfBoundaryZonesXYY == numberOfBoundaryZonesXYZ)
				&& (numberOfBoundaryZonesXYZ == numberOfBoundaryZonesXZX)
				&& (numberOfBoundaryZonesXZX == numberOfBoundaryZonesXZY)
				&& (numberOfBoundaryZonesXZY == numberOfBoundaryZonesXZZ)

				&& (numberOfBoundaryZonesXZZ == numberOfBoundaryZonesYXX)

				&& (numberOfBoundaryZonesYXX == numberOfBoundaryZonesYXY)
				&& (numberOfBoundaryZonesYXY == numberOfBoundaryZonesYXZ)
				&& (numberOfBoundaryZonesYXZ == numberOfBoundaryZonesYYX)
				&& (numberOfBoundaryZonesYYX == numberOfBoundaryZonesYYY)
				&& (numberOfBoundaryZonesYYY == numberOfBoundaryZonesYYZ)
				&& (numberOfBoundaryZonesYYZ == numberOfBoundaryZonesYZX)
				&& (numberOfBoundaryZonesYZX == numberOfBoundaryZonesYZY)
				&& (numberOfBoundaryZonesYZY == numberOfBoundaryZonesYZZ)

				&& (numberOfBoundaryZonesYZZ == numberOfBoundaryZonesZXX)

				&& (numberOfBoundaryZonesZXX == numberOfBoundaryZonesZXY)
				&& (numberOfBoundaryZonesZXY == numberOfBoundaryZonesZXZ)
				&& (numberOfBoundaryZonesZXZ == numberOfBoundaryZonesZYX)
				&& (numberOfBoundaryZonesZYX == numberOfBoundaryZonesZYY)
				&& (numberOfBoundaryZonesZYY == numberOfBoundaryZonesZYZ)
				&& (numberOfBoundaryZonesZYZ == numberOfBoundaryZonesZZX)
				&& (numberOfBoundaryZonesZZX == numberOfBoundaryZonesZZY)
				&& (numberOfBoundaryZonesZZY == numberOfBoundaryZonesZZZ)

		))
			throw exception("Inequal number of tensor3 boundary conditions.",
					errorsEnum::boundaryConditionError);

		if (numberOfBoundaryZonesXXX == 1)
		{
			boundaryConditionValue[0] = tensor3(
					std::stod(boundaryConditionValueStringXXX),
					std::stod(boundaryConditionValueStringXXY),
					std::stod(boundaryConditionValueStringXXZ),
					std::stod(boundaryConditionValueStringXYX),
					std::stod(boundaryConditionValueStringXYY),
					std::stod(boundaryConditionValueStringXYZ),
					std::stod(boundaryConditionValueStringXZX),
					std::stod(boundaryConditionValueStringXZY),
					std::stod(boundaryConditionValueStringXZZ),

					std::stod(boundaryConditionValueStringYXX),
					std::stod(boundaryConditionValueStringYXY),
					std::stod(boundaryConditionValueStringYXZ),
					std::stod(boundaryConditionValueStringYYX),
					std::stod(boundaryConditionValueStringYYY),
					std::stod(boundaryConditionValueStringYYZ),
					std::stod(boundaryConditionValueStringYZX),
					std::stod(boundaryConditionValueStringYZY),
					std::stod(boundaryConditionValueStringYZZ),

					std::stod(boundaryConditionValueStringZXX),
					std::stod(boundaryConditionValueStringZXY),
					std::stod(boundaryConditionValueStringZXZ),
					std::stod(boundaryConditionValueStringZYX),
					std::stod(boundaryConditionValueStringZYY),
					std::stod(boundaryConditionValueStringZYZ),
					std::stod(boundaryConditionValueStringZZX),
					std::stod(boundaryConditionValueStringZZY),
					std::stod(boundaryConditionValueStringZZZ));
		}
		else
		{
			boundaryConditionValue.resize(numberOfBoundaryZonesXXX);

			std::array<std::string, 27> boundaryConditionValueStringV {
					boundaryConditionValueStringXXX,
					boundaryConditionValueStringXXY,
					boundaryConditionValueStringXXZ,
					boundaryConditionValueStringXYX,
					boundaryConditionValueStringXYY,
					boundaryConditionValueStringXYZ,
					boundaryConditionValueStringXZX,
					boundaryConditionValueStringXZY,
					boundaryConditionValueStringXZZ,

					boundaryConditionValueStringYXX,
					boundaryConditionValueStringYXY,
					boundaryConditionValueStringYXZ,
					boundaryConditionValueStringYYX,
					boundaryConditionValueStringYYY,
					boundaryConditionValueStringYYZ,
					boundaryConditionValueStringYZX,
					boundaryConditionValueStringYZY,
					boundaryConditionValueStringYZZ,

					boundaryConditionValueStringZXX,
					boundaryConditionValueStringZXY,
					boundaryConditionValueStringZXZ,
					boundaryConditionValueStringZYX,
					boundaryConditionValueStringZYY,
					boundaryConditionValueStringZYZ,
					boundaryConditionValueStringZZX,
					boundaryConditionValueStringZZY,
					boundaryConditionValueStringZZZ };

			for (std::size_t k = 0; k < tensor3::vsize; ++k)
			{
				auto numBeginIter_k = boundaryConditionValueStringV[k].begin();
				auto numEndIter_k = boundaryConditionValueStringV[k].end();

				for (auto striter = boundaryConditionValueStringV[k].begin();
						striter != boundaryConditionValueStringV[k].end();
						++striter)
				{
					numEndIter_k = striter;

					if (*striter == ';')
						break;
				}

				boundaryConditionValue[0].v_r()[k] = std::stod(
						std::string(numBeginIter_k, numEndIter_k));

				for (std::size_t i = 1; i < numberOfBoundaryZonesXXX; ++i)
				{
					numBeginIter_k = numEndIter_k;
					numBeginIter_k++;

					for (auto striter = numBeginIter_k;
							striter != boundaryConditionValueStringV[k].end();
							++striter)
					{
						numEndIter_k = striter;

						if (*striter == ';')
							break;
					}

					boundaryConditionValue[i].v_r()[k] = std::stod(
							std::string(numBeginIter_k, numEndIter_k));
				}
			}
		}
	}
	else if (boundaryConditionString == "innerSurface")
		throw exception("<<innerSurface>> can not be boundary surface type.",
				errorsEnum::boundaryConditionError);
	else
		throw exception("Unknown type of boundary condition",
				errorsEnum::boundaryConditionError);
}
