/*
 * fabricFunctions.hpp
 *
 *  Created on: 2023/03/16
 *      Author: Maxim Boldyrev
 */

#ifndef FABRICFUNCTIONS_HPP_
#define FABRICFUNCTIONS_HPP_

#include <memory>

#include "chemicalReactionsEnum.hpp"
#include "abstractChemicalKinetics.hpp"
#include "abstractFlowSolver.hpp"
#include "abstractLimiter.hpp"
#include "abstractMatrixSolver.hpp"
#include "abstractTurbulenceGen.hpp"
#include "transportModelEnum.hpp"

namespace schemi
{
std::unique_ptr<abstractLimiter> createLimiter(const std::string & name);

std::pair<std::unique_ptr<abstractMatrixSolver>,
		std::unique_ptr<abstractMatrixSolver>> createMatrixSolver(
		const std::string & name, const std::string & dim,
		const std::size_t iter);

std::unique_ptr<abstractFlowSolver> createFlowSolver(const std::string & name);

std::unique_ptr<abstractTurbulenceGen> createTurbulenceModel(
		const std::string & turbulenceONString,
		const std::string & sourceTypeString);

std::unique_ptr<abstractChemicalKinetics> createChemicalKinetics(
		const homogeneousPhase<cubicCell> & phaseIn,
		const chemicalReactions chemReactFlag, const scalar minimalTimestep);

std::unique_ptr<abstractTransportModel> createTransportModel(
		const std::vector<std::vector<std::string>> & matrixOfSubstancesConditions,
		const scalar constNu, const scalar constD, const scalar constKappa,
		const transportModel model);
}  // namespace schemi

#endif /* FABRICFUNCTIONS_HPP_ */
