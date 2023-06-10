/*
 * abstractMixtureThermodynamics.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractMixtureThermodynamics.hpp"

#include "globalConstants.hpp"

schemi::abstractMixtureThermodynamics::abstractMixtureThermodynamics(
		const scalar Rin, const scalar hPin) noexcept :
		R(Rin), kB(R / NAvogardro), hPlanck(hPin)
{
}

schemi::abstractMixtureThermodynamics::~abstractMixtureThermodynamics() noexcept
{
}
