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
#include "intExpPow.hpp"

std::array<schemi::scalar, 2> schemi::abstractFlowSolver::DavisSignalSpeeds(
		const scalar uL, const scalar aL, const scalar uR,
		const scalar aR) const noexcept
{
	return
	{	std::min(
				uL - aL,
				uR - aR),
		std::max(
				uL + aL,
				uR + aR)};
}

std::array<schemi::scalar, 2> schemi::abstractFlowSolver::EinfeldtSignalSpeeds(
		const scalar rhoL, const scalar uL, const scalar aL, const scalar rhoR,
		const scalar uR, const scalar aR) const noexcept
{
	const auto sRhoL = std::sqrt(rhoL);
	const auto sRhoR = std::sqrt(rhoR);
	const auto sRhoSum = sRhoL + sRhoR;

	const auto uRoe = (uL * sRhoL + uR * sRhoR) / sRhoSum;
	const auto aEnd = std::sqrt(
			(pow<scalar, 2>(aL) * sRhoL + pow<scalar, 2>(aR) * sRhoR) / sRhoSum
					+ 0.5 * sRhoL * sRhoR / pow<scalar, 2>(sRhoSum)
							* pow<scalar, 2>(uR - uL));

	return
	{	std::min(uL - aL, uRoe - aEnd), std::max(uR + aR, uRoe + aEnd)};
}

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
