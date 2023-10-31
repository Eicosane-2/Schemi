/*
 * abstractFlowSolver.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractFlowSolver.hpp"

#include "flowSolverEnum.hpp"
#include "HLLSolver.hpp"
#include "HLLCFSolver.hpp"
#include "HLLCSolver.hpp"
#include "HLLCLMSolver.hpp"
#include "HLLC2pSolver.hpp"
#include "RichtmyerSolver.hpp"
#include "KTSolver.hpp"

schemi::abstractFlowSolver::~abstractFlowSolver() noexcept
{
}

std::unique_ptr<schemi::abstractFlowSolver> schemi::abstractFlowSolver::createFlowSolver(
		const std::string & name, const MPIHandler & par)
{
	flowSolver flowSolwerFlag;
	if (name == "HLLCF")
		flowSolwerFlag = flowSolver::HLLCF;
	else if (name == "HLLC")
		flowSolwerFlag = flowSolver::HLLC;
	else if (name == "HLL")
		flowSolwerFlag = flowSolver::HLL;
	else if (name == "HLLCLM")
		flowSolwerFlag = flowSolver::HLLCLM;
	else if (name == "KT")
		flowSolwerFlag = flowSolver::KT;
	else if (name == "HLLC2p")
		flowSolwerFlag = flowSolver::HLLC2p;
	else if (name == "Richtmyer")
		flowSolwerFlag = flowSolver::Richtmyer;
	else
		throw exception("Unknown flow solver flag.",
				errors::initialisationError);

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
	default:
		return std::make_unique<KTSolver>();
		break;
	}
}
