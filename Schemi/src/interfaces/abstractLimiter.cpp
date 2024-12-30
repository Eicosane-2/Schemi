/*
 * abstractLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "abstractLimiter.hpp"

#include <stdexcept>
#include <map>

#include "exception.hpp"
#include "typeOfTVDLimiterEnum.hpp"
#include "minmodLimiter.hpp"
#include "vanLeerLimiter.hpp"
#include "linearLimiter.hpp"
#include "vanAlbadaLimiter.hpp"
#include "HQUICKLimiter.hpp"
#include "vanLeer2Limiter.hpp"
#include "superbeeLimiter.hpp"
#include "vanAlbada2Limiter.hpp"
#include "minmod2Limiter.hpp"
#include "SewbyLimiter.hpp"
#include "zeroLimiter.hpp"

schemi::abstractLimiter::~abstractLimiter() noexcept
{
}

std::unique_ptr<schemi::abstractLimiter> schemi::abstractLimiter::createLimiter(
		const std::string & name)
{
	std::map<std::string, typeOfTVDLimiter> limiterType;
	limiterType.insert( { "minmod", typeOfTVDLimiter::minmod });
	limiterType.insert( { "vanLeer", typeOfTVDLimiter::vanLeer });
	limiterType.insert( { "linear", typeOfTVDLimiter::linear });
	limiterType.insert( { "zero", typeOfTVDLimiter::zero });
	limiterType.insert( { "vanAlbada", typeOfTVDLimiter::vanAlbada });
	limiterType.insert( { "HQUICK", typeOfTVDLimiter::HQUICK });
	limiterType.insert( { "vanLeer2", typeOfTVDLimiter::vanLeer2 });
	limiterType.insert( { "superbee", typeOfTVDLimiter::superbee });
	limiterType.insert( { "vanAlbada2", typeOfTVDLimiter::vanAlbada2 });
	limiterType.insert( { "minmod2", typeOfTVDLimiter::minmod2 });
	limiterType.insert( { "Sewby", typeOfTVDLimiter::Sewby });

	typeOfTVDLimiter limiterFlag;
	try
	{
		limiterFlag = limiterType.at(name);
	} catch (const std::out_of_range&)
	{
		throw exception("Unknown TVD limiter flag.",
				errors::initialisationError);
	}

	switch (limiterFlag)
	{
	case typeOfTVDLimiter::minmod:
		return std::make_unique<minmodLimiter>();
		break;
	case typeOfTVDLimiter::vanLeer:
		return std::make_unique<vanLeerLimiter>();
		break;
	case typeOfTVDLimiter::linear:
		return std::make_unique<linearLimiter>();
		break;
	case typeOfTVDLimiter::vanAlbada:
		return std::make_unique<vanAlbadaLimiter>();
		break;
	case typeOfTVDLimiter::HQUICK:
		return std::make_unique<HQUICKLimiter>();
		break;
	case typeOfTVDLimiter::vanLeer2:
		return std::make_unique<vanLeer2Limiter>();
		break;
	case typeOfTVDLimiter::superbee:
		return std::make_unique<superbeeLimiter>();
		break;
	case typeOfTVDLimiter::vanAlbada2:
		return std::make_unique<vanAlbada2Limiter>();
		break;
	case typeOfTVDLimiter::minmod2:
		return std::make_unique<minmod2Limiter>();
		break;
	case typeOfTVDLimiter::Sewby:
		return std::make_unique<SewbyLimiter>();
		break;
	[[unlikely]] default:
		return std::make_unique<zeroLimiter>();
		break;
	}
}
