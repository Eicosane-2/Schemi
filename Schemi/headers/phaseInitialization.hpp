/*
 * phaseInitialization.hpp
 *
 *  Created on: 2022/10/19
 *      Author: Maxim Boldyrev
 *
 *      Initialization of phase data.
 */

#ifndef PHASEINITIALIZATION_HPP_
#define PHASEINITIALIZATION_HPP_

#include <tuple>
#include <memory>

#include "enthalpyFlowEnum.hpp"
#include "homogeneousPhase.hpp"
#include "MPIHandler.hpp"

namespace schemi
{
std::tuple<std::unique_ptr<homogeneousPhase<cubicCell>>, enthalpyFlowEnum, bool> phaseInitialization(
		std::size_t numberOfComponents, std::size_t numberOfZones,
		const mesh & meshReference,
		const std::vector<boundaryConditionType> & commonConditions,
		const MPIHandler & parallelism, const std::string & turbulenceONString,
		const std::string & sourceTypeString,
		const std::string & universalGasConstant,
		const std::string & equationOfState, const std::size_t readDataPoint);
}  // namespace schemi

#endif /* PHASEINITIALIZATION_HPP_ */
