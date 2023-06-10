/*
 * conservativeFlows.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "conservativeFlows.hpp"

schemi::conservativeFlows::conservativeFlows(const mesh & meshRef,
		const std::size_t numberOfcomponents) noexcept :
		density { numberOfcomponents + 1, surfaceField<vector>(meshRef,
				vector(0)) },

		momentum { meshRef, tensor(0) },

		totalEnergy { meshRef, vector(0) },

		rhokTurb { meshRef, vector { 0 } },

		rhoepsTurb { meshRef, vector(0) },

		rhoaTurb { meshRef, tensor(0) },

		rhobTurb { meshRef, vector(0) }
{
}
