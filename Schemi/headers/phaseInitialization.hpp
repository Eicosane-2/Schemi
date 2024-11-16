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
#include "MPIHandler.hpp"
#include "homogeneousPhase.hpp"
#include "zone.hpp"

namespace schemi
{
std::tuple<std::unique_ptr<homogeneousPhase<cubicCell>>, enthalpyFlow, bool> phaseInitialization(
		std::size_t numberOfComponents,
		const std::vector<std::unique_ptr<zone>> & numberOfZones,
		const mesh & meshReference,
		const std::vector<boundaryConditionType> & commonConditions,
		const MPIHandler & parallelism, const std::string & turbulenceONString,
		const std::string & sourceTypeString,
		const std::string & universalGasConstant,
		const std::string & equationOfState,
		const std::pair<std::size_t, std::string> & readDataPoint);
}  // namespace schemi

#endif /* PHASEINITIALIZATION_HPP_ */
