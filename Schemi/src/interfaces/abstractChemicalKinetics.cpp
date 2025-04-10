/*
 * abstractChemicalKinetics.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "abstractChemicalKinetics.hpp"

#include "chemicalKineticsChlorumDissociation.hpp"
#include "chemicalKineticsChlorumHydrogeniumDissociation.hpp"
#include "chemicalKineticsH2Cl2Combustion.hpp"
#include "chemicalKineticsH2O2Combustion.hpp"
#include "chemicalKineticsNO2Disproportionation.hpp"
#include "chemicalKineticsNoReaction.hpp"
#include "chemicalKineticsRober.hpp"

schemi::chemicalKinetics::abstractChemicalKinetics::abstractChemicalKinetics(
		const bool flag, const scalar mt) noexcept :
		minTimestep(mt), solvers(), chemicalReaction(flag)
{
	solvers.insert( { "Gauss-Seidel", iterativeSolver::GaussSeidel });
	solvers.insert(
			{ "conjugateGradient", iterativeSolver::ConjugateGradient });
	solvers.insert( { "JacobiConjugateGradient",
			iterativeSolver::JacobiConjugateGradient });
	solvers.insert( { "Jacobi", iterativeSolver::Jacobi });
	solvers.insert( { "GaussElimination", iterativeSolver::GaussElimination });
	solvers.insert( { "no", iterativeSolver::noSolver });
}

schemi::chemicalKinetics::abstractChemicalKinetics::~abstractChemicalKinetics() noexcept
{
}

std::unique_ptr<schemi::chemicalKinetics::abstractChemicalKinetics> schemi::chemicalKinetics::abstractChemicalKinetics::createChemicalKinetics(
		const homogeneousPhase<cubicCell> & phaseIn,
		const chemicalReactions chemReactFlag,
		const scalar minimalTimestep) noexcept
{
	switch (chemReactFlag)
	{
	case chemicalReactions::Cl2Dissociation:
		return std::make_unique<chemicalKinetics::ChlorumDissociation>(phaseIn,
				minimalTimestep);
		break;
	case chemicalReactions::Cl2H2Dissociation:
		return std::make_unique<chemicalKinetics::ChlorumHydrogeniumDissociation>(
				phaseIn, minimalTimestep);
		break;
	case chemicalReactions::H2Cl2Combustion:
		return std::make_unique<chemicalKinetics::H2Cl2Combustion>(phaseIn,
				minimalTimestep);
		break;
	case chemicalReactions::NO2Disproportionation:
		return std::make_unique<chemicalKinetics::NO2Disproportionation>(
				phaseIn, minimalTimestep);
		break;
	case chemicalReactions::H2O2Combustion:
		return std::make_unique<chemicalKinetics::H2O2Combustion>(phaseIn,
				minimalTimestep);
		break;
	case chemicalReactions::Rober:
		return std::make_unique<chemicalKinetics::Rober>(phaseIn,
				minimalTimestep);
		break;
	case chemicalReactions::noReaction:
	default:
		return std::make_unique<chemicalKinetics::NoReaction>();
		break;
	}
}

void schemi::chemicalKinetics::abstractChemicalKinetics::normalize(
		std::valarray<scalar> & res) noexcept
{
	for (auto & i : res)
		if (std::abs(i) < std::numeric_limits<scalar>::epsilon())
			i = 0;
}
