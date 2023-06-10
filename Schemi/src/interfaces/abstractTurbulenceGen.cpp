/*
 * abstractTurbulenceGen.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractTurbulenceGen.hpp"

schemi::abstractTurbulenceGen::abstractTurbulenceGen(const bool turb_in,
		const turbulenceModelEnum tm_in) noexcept :
		turbulence(turb_in), model(tm_in)
{
}

schemi::abstractTurbulenceGen::~abstractTurbulenceGen() noexcept
{
}
