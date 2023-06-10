/*
 * conservativeFlows.hpp
 *
 *  Created on: 2019/12/22
 *      Author: Maxim Boldyrev
 *
 *      Structure for storing conservative fields flows.
 */

#ifndef CONSERVATIVEFLOWS_HPP_
#define CONSERVATIVEFLOWS_HPP_

#include <vector>

#include "mesh.hpp"
#include "tensor.hpp"
#include "vector.hpp"
#include "surfaceField.hpp"

namespace schemi
{
struct conservativeFlows
{
	std::vector<surfaceField<vector>> density;
	surfaceField<tensor> momentum;
	surfaceField<vector> totalEnergy;
	surfaceField<vector> rhokTurb;
	surfaceField<vector> rhoepsTurb;
	surfaceField<tensor> rhoaTurb;
	surfaceField<vector> rhobTurb;

	conservativeFlows(const mesh & meshRef,
			const std::size_t numberOfcomponents) noexcept;
};
}  // namespace schemi

#endif /* CONSERVATIVEFLOWS_HPP_ */
