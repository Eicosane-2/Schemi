/*
 * abstractChemicalKinetics.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractChemicalKinetics.hpp"

schemi::abstractChemicalKinetics::abstractChemicalKinetics(
		const bool flag) noexcept :
		chemicalReaction(flag)
{
}

schemi::abstractChemicalKinetics::~abstractChemicalKinetics() noexcept
{
}
