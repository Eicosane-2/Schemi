/*
 * abstractLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "abstractLimiter.hpp"

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
		const std::string_view name)
{
	typeOfTVDLimiter limiterFlag;
	if (name == "minmod")
		limiterFlag = typeOfTVDLimiter::minmod;
	else if (name == "vanLeer")
		limiterFlag = typeOfTVDLimiter::vanLeer;
	else if (name == "linear")
		limiterFlag = typeOfTVDLimiter::linear;
	else if (name == "zero")
		limiterFlag = typeOfTVDLimiter::zero;
	else if (name == "vanAlbada")
		limiterFlag = typeOfTVDLimiter::vanAlbada;
	else if (name == "HQUICK")
		limiterFlag = typeOfTVDLimiter::HQUICK;
	else if (name == "vanLeer2")
		limiterFlag = typeOfTVDLimiter::vanLeer2;
	else if (name == "superbee")
		limiterFlag = typeOfTVDLimiter::superbee;
	else if (name == "vanAlbada2")
		limiterFlag = typeOfTVDLimiter::vanAlbada2;
	else if (name == "minmod2")
		limiterFlag = typeOfTVDLimiter::minmod2;
	else if (name == "Sewby")
		limiterFlag = typeOfTVDLimiter::Sewby;
	else
		throw exception("Unknown TVD limiter flag.",
				errors::initialisationError);

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
	default:
		return std::make_unique<zeroLimiter>();
		break;
	}
}
