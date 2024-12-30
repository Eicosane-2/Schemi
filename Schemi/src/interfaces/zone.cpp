/*
 * zone.cpp
 *
 *  Created on: 2024/10/09
 *      Author: Maxim Boldyrev
 */

#include "zone.hpp"

#include <iostream>
#include <fstream>

#include "exception.hpp"
#include "outerZone.hpp"
#include "zoneInsideCylinder.hpp"
#include "zoneInsideSphere.hpp"
#include "zoneUnderPlane.hpp"
#include "zoneUnderPlanes.hpp"
#include "globalConstants.hpp"

schemi::zone::~zone() noexcept
{
}

std::vector<std::unique_ptr<schemi::zone>> schemi::zone::zonesArray()
{
	std::ifstream zoneDescr("./set/zonesDescription.txt");
	if (zoneDescr.is_open())
		std::cout << "./set/zonesDescription.txt is opened." << std::endl;
	else
		[[unlikely]]
		throw std::ifstream::failure("./set/zonesDescription.txt not found.");

	std::string skipBuffer, outZone;

	std::size_t nzones;

	zoneDescr >> skipBuffer >> nzones;
	zoneDescr >> skipBuffer >> outZone;

	bool IsOuterZone;
	try
	{
		IsOuterZone = onOffMap.at(outZone);
	} catch (const std::out_of_range&)
	{
		throw exception("Unknown flag for outer zone description.",
				errors::initialisationError);
	}

	zoneDescr >> skipBuffer;

	std::vector<std::unique_ptr<zone>> zonesCollection(nzones);

	for (std::size_t z = 0; z < zonesCollection.size(); ++z)
	{
		std::string zoneType;

		zoneDescr >> zoneType;

		if (zoneType == "plane")
		{
			scalar nX, nY, nZ, pX, pY, pZ;
			zoneDescr >> nX >> nY >> nZ;
			zoneDescr >> skipBuffer;
			zoneDescr >> pX >> pY >> pZ;

			zonesCollection[z] = std::make_unique<zoneUnderPlane>(
					vector(nX, nY, nZ), vector(pX, pY, pZ));
		}
		else if (zoneType == "planes")
		{
			std::size_t nPlanes;
			zoneDescr >> nPlanes;

			std::vector<vector> planesN, planesP;

			for (std::size_t pl = 0; pl < nPlanes; ++pl)
			{
				scalar nX, nY, nZ, pX, pY, pZ;
				zoneDescr >> nX >> nY >> nZ;
				zoneDescr >> skipBuffer;
				zoneDescr >> pX >> pY >> pZ;

				planesN.push_back(vector(nX, nY, nZ));
				planesP.push_back(vector(pX, pY, pZ));
			}

			zonesCollection[z] = std::make_unique<zoneUnderPlanes>(planesN,
					planesP);
		}
		else if (zoneType == "sphere")
		{
			scalar pX, pY, pZ, r;

			zoneDescr >> pX >> pY >> pZ;
			zoneDescr >> r;

			zonesCollection[z] = std::make_unique<zoneInsideSphere>(
					vector(pX, pY, pZ), r);
		}
		else if (zoneType == "cylinder")
		{
			scalar nX, nY, nZ, pX1, pY1, pZ1, pX2, pY2, pZ2, r;

			zoneDescr >> nX >> nY >> nZ;
			zoneDescr >> skipBuffer;
			zoneDescr >> pX1 >> pY1 >> pZ1;
			zoneDescr >> skipBuffer;
			zoneDescr >> pX2 >> pY2 >> pZ2;
			zoneDescr >> r;

			zonesCollection[z] = std::make_unique<zoneInsideCylinder>(
					vector(nX, nY, nZ), vector(pX1, pY1, pZ1),
					vector(pX2, pY2, pZ2), r);
		}
		else
			[[unlikely]]
			throw exception("Unknown zone type.", errors::initialisationError);
	}

	if (IsOuterZone)
		zonesCollection.push_back(std::make_unique<outerZone>());

	return zonesCollection;
}
