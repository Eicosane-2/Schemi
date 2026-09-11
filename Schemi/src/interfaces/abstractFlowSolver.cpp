/*
 * abstractFlowSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractFlowSolver.hpp"

#include <map>

#include "flowSolverEnum.hpp"
#include "HLLSolver.hpp"
#include "HLLCFSolver.hpp"
#include "HLLCSolver.hpp"
#include "HLLCLMSolver.hpp"
#include "HLLC2pSolver.hpp"
#include "RichtmyerSolver.hpp"
#include "KTSolver.hpp"
#include "HLLCKTSolver.hpp"

schemi::abstractFlowSolver::~abstractFlowSolver() noexcept
{
}

std::unique_ptr<schemi::abstractFlowSolver> schemi::abstractFlowSolver::createFlowSolver(
		const std::string & name, const MPIHandler & par)
{
	std::map<std::string, flowSolver> flowSolverType;
	flowSolverType.insert( { "HLLCF", flowSolver::HLLCF });
	flowSolverType.insert( { "HLLC", flowSolver::HLLC });
	flowSolverType.insert( { "HLL", flowSolver::HLL });
	flowSolverType.insert( { "HLLCLM", flowSolver::HLLCLM });
	flowSolverType.insert( { "KT", flowSolver::KT });
	flowSolverType.insert( { "HLLC2p", flowSolver::HLLC2p });
	flowSolverType.insert( { "Richtmyer", flowSolver::Richtmyer });
	flowSolverType.insert( { "HLLCKT", flowSolver::HLLCKT });

	flowSolver flowSolwerFlag;
	try
	{
		flowSolwerFlag = flowSolverType.at(name);
	} catch (const std::out_of_range&)
	{
		throw exception("Unknown flow solver flag.",
				errors::initialisationError);
	}

	switch (flowSolwerFlag)
	{
	case flowSolver::HLL:
		return std::make_unique<HLLSolver>();
		break;
	case flowSolver::HLLCF:
		return std::make_unique<HLLCFSolver>();
		break;
	case flowSolver::HLLC:
		return std::make_unique<HLLCSolver>();
		break;
	case flowSolver::HLLCLM:
		return std::make_unique<HLLCLMSolver>();
		break;
	case flowSolver::HLLC2p:
		return std::make_unique<HLLC2pSolver>();
		break;
	case flowSolver::Richtmyer:
		return std::make_unique<RichtmyerSolver>(par);
		break;
	case flowSolver::HLLCKT:
		return std::make_unique<HLLCKTSolver>();
		break;
	[[unlikely]] default:
		return std::make_unique<KTSolver>();
		break;
	}
}
