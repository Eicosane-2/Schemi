/*
 * chemicalKineticsNoReaction.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsNoReaction.hpp"

schemi::chemicalKinetics::NoReaction::NoReaction() noexcept :
		abstractChemicalKinetics(false, 0.0)
{
}

void schemi::chemicalKinetics::NoReaction::solveChemicalKinetics(
		homogeneousPhase<cubicCell>&) const
{
}
