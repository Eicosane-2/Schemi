/*
 * chemicalKineticsNoReaction.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "chemicalKineticsNoReaction.hpp"

schemi::chemicalKineticsNoReaction::chemicalKineticsNoReaction() noexcept :
		abstractChemicalKinetics(false)
{
}

void schemi::chemicalKineticsNoReaction::solveChemicalKinetics(
		homogeneousPhase<cubicCell>&) const noexcept
{
}
