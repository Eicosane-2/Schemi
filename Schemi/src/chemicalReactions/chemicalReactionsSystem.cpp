/*
 * chemicalReactionsSystem.cpp
 *
 *  Created on: 2026/06/30
 *      Author: Maxim Boldyrev
 */

#include "chemicalReactionsSystem.hpp"

#include <iostream>
#include <fstream>
#include <string>

#include "exception.hpp"
#include "globalConstants.hpp"
#include "intExpPow.hpp"

std::valarray<schemi::scalar> schemi::chemicalKinetics::chemicalReactionsSystem::GaussElemination(
		std::valarray<scalar> & A, std::valarray<scalar> & b) const
{
	const auto N { b.size() };

	std::valarray<scalar> w(N);

	for (std::size_t k = 0; k < N - 1; ++k)
	{
		for (std::size_t i = k + 1; i < N; ++i)
		{
			const auto ik = i * N + k;
			const auto kk = k * N + k;
			const auto Ratio = A[ik] / (A[kk] + stabilizator);
			for (std::size_t j = k + 1; j < N; ++j)
			{
				const auto ij = i * N + j;
				const auto kj = k * N + j;
				A[ij] -= Ratio * A[kj];
			}
			b[i] -= Ratio * b[k];
		}
	}
	w[N - 1] = b[N - 1] / (A[pow<scalar, 2>(N) - 1] + stabilizator);
	for (std::size_t i = N - 2;; --i)
	{
		scalar Term = 0.;
		for (std::size_t j = i + 1; j < N; ++j)
		{
			const auto ij = i * N + j;
			Term += A[ij] * w[j];
		}
		const auto ii = i * N + i;
		w[i] = (b[i] - Term) / (A[ii] + stabilizator);

		if (i == 0)
			break;
	}

	return w;
}

void schemi::chemicalKinetics::chemicalReactionsSystem::renormalization(
		const scalar sumMassFracOld,
		std::valarray<scalar> & massFractions) const
{
	for (const auto & mfr : massFractions)
		if (mfr < 0)
			throw exception("New mass fraction is negative.",
					errors::positivnessError);

	const auto sumMassFracNew = massFractions.sum();
	if (std::abs(sumMassFracNew - sumMassFracOld) > massFracTolerance)
		throw exception("Mass fraction difference too big.",
				errors::positivnessError);
	else
	{
		massFractions /= sumMassFracNew;
		massFractions *= sumMassFracOld;
	}
}

schemi::chemicalKinetics::chemicalReactionsSystem::chemicalReactionsSystem(
		const abstractMixtureThermodynamics & thermIn) :
		therm(thermIn)
{
	std::string buffer, swt;

	std::ifstream chemKinDescription { "./set/chemicalReactionsSystem.txt" };
	if (chemKinDescription.is_open())
		std::cout << "./set/chemicalReactionsSystem.txt is opened."
				<< std::endl;
	else
	{
		[[unlikely]]
		throw std::ifstream::failure(
				"./set/chemicalReactionsSystem.txt not found.");
	}

	chemKinDescription >> buffer >> swt;

	chemicalReactions = onOffMap.at(swt);

	if (chemicalReactions)
	{
		std::string chemReactionDescription;

		chemKinDescription >> buffer;

		chemKinDescription.ignore();

		while (std::getline(chemKinDescription, chemReactionDescription))
			reactionsParameters.emplace_back(
					irreversibleReaction(chemReactionDescription));

		std::vector<std::string> reactingComponents;

		for (const auto & reactionsParameters_i : reactionsParameters)
		{
			for (const auto & reagent : reactionsParameters_i.reagentSubstNamesGet())
			{
				if (reagent == "M0")
					continue;

				bool found(false);
				for (const auto & component : reactingComponents)
					if (component == reagent)
					{
						found = true;
						break;
					}

				if (!found)
					reactingComponents.push_back(reagent);
			}

			for (const auto & product : reactionsParameters_i.productSubstNamesGet())
			{
				if (product == "M0")
					continue;

				bool found(false);
				for (const auto & component : reactingComponents)
					if (component == product)
					{
						found = true;
						break;
					}

				if (!found)
					reactingComponents.push_back(product);
			}
		}

		const auto & names = therm.getSubstancesNames();

		for (const auto & reactComp : reactingComponents)
		{
			bool found(false);
			for (std::size_t k = 0; k < names.size(); ++k)
				if (names[k] == reactComp)
				{
					reactingComponentsMatching[reactComp] = k;
					found = true;
					break;
				}

			if (!found)
				throw exception("Substance name is not found in thermodynamic.",
						errors::initialisationError);
		}

		std::cout << "Reagents: ";
		for (const auto & comp : reactingComponentsMatching)
		{
			reactingComponentsIndexes.push_back(comp.second);
			std::cout << comp.first << ' ' << comp.second + 1 << ", ";
		}
		std::cout << '.' << std::endl;

		std::sort(reactingComponentsIndexes.begin(),
				reactingComponentsIndexes.end());

		matrixPrototype.resize(
				pow<std::size_t, 2>(reactingComponentsMatching.size()));

		for (std::size_t iCell = 0; iCell < matrixPrototype.size(); ++iCell)
		{
			const auto matrInd_i = iCell / reactingComponentsMatching.size();
			const auto matrInd_j = iCell
					- reactingComponentsMatching.size() * matrInd_i;
			const auto nR = reactingComponentsIndexes[matrInd_i];
			const auto nC = reactingComponentsIndexes[matrInd_j];

			matrixPrototype[iCell].indexes = { nR + 1, nC + 1 };

			if (matrInd_i == matrInd_j)
			{
				matrixPrototype[iCell].nullCell = false;
				matrixPrototype[iCell].diagonalCell = true;

				matrixPrototype[iCell].molMass.first = therm.Mv()[nR];
			}
			else
			{
				matrixPrototype[iCell].nullCell = true;
				matrixPrototype[iCell].diagonalCell = false;
			}

			const auto & compNameR = therm.getSubstancesNames()[nR];
			const auto & compNameC = therm.getSubstancesNames()[nC];

			for (const auto & r : reactionsParameters)
			{
				const bool cRow = std::find(r.reagentSubstNamesGet().cbegin(),
						r.reagentSubstNamesGet().cend(), compNameR)
						!= r.reagentSubstNamesGet().cend();
				const bool cCol = std::find(r.reagentSubstNamesGet().cbegin(),
						r.reagentSubstNamesGet().cend(), compNameC)
						!= r.reagentSubstNamesGet().cend();

				if (cRow && cCol)
				{
					matrixPrototype[iCell].nullCell = false;

					matrixPrototype[iCell].reactParams.push_back(
							r.reactionParametersGet());

					std::size_t numberOfIndivComps =
							r.reagentSubstNamesGet().size();
					if (std::find(r.reagentSubstNamesGet().cbegin(),
							r.reagentSubstNamesGet().cend(), "M0")
							!= r.reagentSubstNamesGet().cend())
						numberOfIndivComps--;

					matrixPrototype[iCell].reactWeight.push_back(
							1. / scalar(numberOfIndivComps));

					std::vector<std::pair<std::size_t, int>> compIndexes;
					for (std::size_t jComp = 0;
							jComp < r.reagentSubstNamesGet().size(); ++jComp)
					{
						const auto & rComp = r.reagentSubstNamesGet()[jComp];
						if (rComp == "M0")
							compIndexes.emplace_back(
									std::make_pair(std::size_t(0),
											r.reagentCoeffsGet()[jComp]));
						else if (rComp == compNameC)
						{
							const auto compC = r.reagentCoeffsGet()[jComp];
							if (compC != 1)
								compIndexes.emplace_back(
										std::make_pair(
												reactingComponentsMatching.at(
														rComp) + 1, compC - 1));
						}
						else
							compIndexes.emplace_back(
									std::make_pair(
											reactingComponentsMatching.at(rComp)
													+ 1,
											r.reagentCoeffsGet()[jComp]));
					}
					matrixPrototype[iCell].comp.push_back(compIndexes);

					matrixPrototype[iCell].molMass.second.push_back(
							therm.Mv()[nC]);

					for (std::size_t jComp = 0;
							jComp < r.reagentSubstNamesGet().size(); ++jComp)
						if (r.reagentSubstNamesGet()[jComp] == compNameR)
						{
							matrixPrototype[iCell].coeffs.push_back(
									r.reagentCoeffsGet()[jComp]);
							break;
						}
				}
			}

			for (const auto & r : reactionsParameters)
			{
				const bool cRow = std::find(r.productSubstNamesGet().cbegin(),
						r.productSubstNamesGet().cend(), compNameR)
						!= r.productSubstNamesGet().cend();
				const bool cCol = std::find(r.reagentSubstNamesGet().cbegin(),
						r.reagentSubstNamesGet().cend(), compNameC)
						!= r.reagentSubstNamesGet().cend();

				if (cRow && cCol)
				{
					matrixPrototype[iCell].nullCell = false;

					matrixPrototype[iCell].reactParams.push_back(
							r.reactionParametersGet());

					std::size_t numberOfIndivComps =
							r.reagentSubstNamesGet().size();
					if (std::find(r.reagentSubstNamesGet().cbegin(),
							r.reagentSubstNamesGet().cend(), "M0")
							!= r.reagentSubstNamesGet().cend())
						numberOfIndivComps--;

					matrixPrototype[iCell].reactWeight.push_back(
							1. / scalar(numberOfIndivComps));

					std::vector<std::pair<std::size_t, int>> compIndexes;
					for (std::size_t jComp = 0;
							jComp < r.reagentSubstNamesGet().size(); ++jComp)
					{
						const auto & rComp = r.reagentSubstNamesGet()[jComp];
						if (rComp == "M0")
							compIndexes.emplace_back(
									std::make_pair(std::size_t(0),
											r.reagentCoeffsGet()[jComp]));
						else if (rComp == compNameC)
						{
							const auto compC = r.reagentCoeffsGet()[jComp];
							if (compC != 1)
								compIndexes.emplace_back(
										std::make_pair(
												reactingComponentsMatching.at(
														rComp) + 1, compC - 1));
						}
						else
							compIndexes.emplace_back(
									std::make_pair(
											reactingComponentsMatching.at(rComp)
													+ 1,
											r.reagentCoeffsGet()[jComp]));
					}
					matrixPrototype[iCell].comp.push_back(compIndexes);

					matrixPrototype[iCell].molMass.second.push_back(
							therm.Mv()[nC]);

					for (std::size_t jComp = 0;
							jComp < r.productSubstNamesGet().size(); ++jComp)
						if (r.productSubstNamesGet()[jComp] == compNameR)
						{
							matrixPrototype[iCell].coeffs.push_back(
									-r.productCoeffsGet()[jComp]);
							break;
						}
				}
			}

			/*if (matrInd_i == matrInd_j)
			 {
			 matrixPrototype[iCell].nullCell = false;
			 matrixPrototype[iCell].diagonalCell = true;

			 const auto & compName = therm.getSubstancesNames()[nC];

			 matrixPrototype[iCell].molMass.first = therm.Mv()[nC];

			 for (const auto & r : reactionsParameters)
			 {
			 if (std::find(r.reagentSubstNamesGet().cbegin(),
			 r.reagentSubstNamesGet().cend(), compName)
			 != r.reagentSubstNamesGet().cend())
			 {
			 matrixPrototype[iCell].reactParams.push_back(
			 r.reactionParametersGet());

			 std::size_t numberOfIndivComps =
			 r.reagentSubstNamesGet().size();
			 if (std::find(r.reagentSubstNamesGet().cbegin(),
			 r.reagentSubstNamesGet().cend(), "M0")
			 != r.reagentSubstNamesGet().cend())
			 numberOfIndivComps--;

			 matrixPrototype[iCell].reactWeight.push_back(
			 1. / scalar(numberOfIndivComps));

			 std::vector<std::pair<std::size_t, int>> compIndexes;
			 for (std::size_t jComp = 0;
			 jComp < r.reagentSubstNamesGet().size();
			 ++jComp)
			 {
			 const auto & rComp = r.reagentSubstNamesGet()[jComp];
			 if (rComp == "M0")
			 compIndexes.emplace_back(
			 std::make_pair(std::size_t(0),
			 r.reagentCoeffsGet()[jComp]));
			 else if (rComp == compName)
			 {
			 const auto compC = r.reagentCoeffsGet()[jComp];
			 if (compC != 1)
			 compIndexes.emplace_back(
			 std::make_pair(
			 reactingComponentsMatching.at(
			 rComp) + 1,
			 compC - 1));
			 }
			 else
			 compIndexes.emplace_back(
			 std::make_pair(
			 reactingComponentsMatching.at(
			 rComp) + 1,
			 r.reagentCoeffsGet()[jComp]));
			 }
			 matrixPrototype[iCell].comp.push_back(compIndexes);

			 matrixPrototype[iCell].molMass.second.push_back(
			 therm.Mv()[nR]);

			 for (std::size_t jComp = 0;
			 jComp < r.reagentSubstNamesGet().size();
			 ++jComp)
			 if (r.reagentSubstNamesGet()[jComp] == compName)
			 {
			 matrixPrototype[iCell].coeffs.push_back(
			 r.reagentCoeffsGet()[jComp]);
			 break;
			 }
			 }
			 }

			 for (const auto & r : reactionsParameters)
			 {
			 if (std::find(r.productSubstNamesGet().cbegin(),
			 r.productSubstNamesGet().cend(), compName)
			 != r.productSubstNamesGet().cend())
			 {
			 matrixPrototype[iCell].reactParams.push_back(
			 r.reactionParametersGet());

			 std::size_t numberOfIndivComps =
			 r.reagentSubstNamesGet().size();
			 if (std::find(r.reagentSubstNamesGet().cbegin(),
			 r.reagentSubstNamesGet().cend(), "M0")
			 != r.reagentSubstNamesGet().cend())
			 numberOfIndivComps--;

			 matrixPrototype[iCell].reactWeight.push_back(
			 1. / scalar(numberOfIndivComps));

			 std::vector<std::pair<std::size_t, int>> compIndexes;
			 for (std::size_t jComp = 0;
			 jComp < r.reagentSubstNamesGet().size();
			 ++jComp)
			 {
			 const auto & rComp = r.reagentSubstNamesGet()[jComp];
			 if (rComp == "M0")
			 compIndexes.emplace_back(
			 std::make_pair(std::size_t(0),
			 r.reagentCoeffsGet()[jComp]));
			 else if (rComp == compName)
			 {
			 const auto compC = r.reagentCoeffsGet()[jComp];
			 if (compC != 1)
			 compIndexes.emplace_back(
			 std::make_pair(
			 reactingComponentsMatching.at(
			 rComp) + 1,
			 compC - 1));
			 }
			 else
			 compIndexes.emplace_back(
			 std::make_pair(
			 reactingComponentsMatching.at(
			 rComp) + 1,
			 r.reagentCoeffsGet()[jComp]));
			 }
			 matrixPrototype[iCell].comp.push_back(compIndexes);

			 matrixPrototype[iCell].molMass.second.push_back(
			 therm.Mv()[nR]);

			 for (std::size_t jComp = 0;
			 jComp < r.reagentSubstNamesGet().size();
			 ++jComp)
			 if (r.reagentSubstNamesGet()[jComp] == compName)
			 {
			 matrixPrototype[iCell].coeffs.push_back(
			 -r.reagentCoeffsGet()[jComp]);
			 break;
			 }
			 }
			 }
			 }
			 else
			 {
			 matrixPrototype[iCell].nullCell = true;
			 matrixPrototype[iCell].diagonalCell = false;

			 const auto & compNameR = therm.getSubstancesNames()[nR];
			 const auto & compNameC = therm.getSubstancesNames()[nC];

			 for (const auto & r : reactionsParameters)
			 {
			 const bool cRow = std::find(
			 r.reagentSubstNamesGet().cbegin(),
			 r.reagentSubstNamesGet().cend(), compNameR)
			 != r.reagentSubstNamesGet().cend();
			 const bool cCol = std::find(
			 r.reagentSubstNamesGet().cbegin(),
			 r.reagentSubstNamesGet().cend(), compNameC)
			 != r.reagentSubstNamesGet().cend();

			 if (cRow && cCol)
			 {
			 matrixPrototype[iCell].nullCell = false;

			 matrixPrototype[iCell].reactParams.push_back(
			 r.reactionParametersGet());

			 std::size_t numberOfIndivComps =
			 r.reagentSubstNamesGet().size();
			 if (std::find(r.reagentSubstNamesGet().cbegin(),
			 r.reagentSubstNamesGet().cend(), "M0")
			 != r.reagentSubstNamesGet().cend())
			 numberOfIndivComps--;

			 matrixPrototype[iCell].reactWeight.push_back(
			 1. / scalar(numberOfIndivComps));

			 std::vector<std::pair<std::size_t, int>> compIndexes;
			 for (std::size_t jComp = 0;
			 jComp < r.reagentSubstNamesGet().size();
			 ++jComp)
			 {
			 const auto & rComp = r.reagentSubstNamesGet()[jComp];
			 if (rComp == "M0")
			 compIndexes.emplace_back(
			 std::make_pair(std::size_t(0),
			 r.reagentCoeffsGet()[jComp]));
			 else if (rComp == compNameC)
			 {
			 const auto compC = r.reagentCoeffsGet()[jComp];
			 if (compC != 1)
			 compIndexes.emplace_back(
			 std::make_pair(
			 reactingComponentsMatching.at(
			 rComp) + 1,
			 compC - 1));
			 }
			 else
			 compIndexes.emplace_back(
			 std::make_pair(
			 reactingComponentsMatching.at(
			 rComp) + 1,
			 r.reagentCoeffsGet()[jComp]));
			 }
			 matrixPrototype[iCell].comp.push_back(compIndexes);

			 matrixPrototype[iCell].molMass.second.push_back(
			 therm.Mv()[nR]);

			 for (std::size_t jComp = 0;
			 jComp < r.reagentSubstNamesGet().size();
			 ++jComp)
			 if (r.reagentSubstNamesGet()[jComp] == compNameR)
			 {
			 matrixPrototype[iCell].coeffs.push_back(
			 r.reagentCoeffsGet()[jComp]);
			 break;
			 }
			 }
			 }

			 for (const auto & r : reactionsParameters)
			 {
			 const bool cRow = std::find(
			 r.productSubstNamesGet().cbegin(),
			 r.productSubstNamesGet().cend(), compNameR)
			 != r.productSubstNamesGet().cend();
			 const bool cCol = std::find(
			 r.reagentSubstNamesGet().cbegin(),
			 r.reagentSubstNamesGet().cend(), compNameC)
			 != r.reagentSubstNamesGet().cend();

			 if (cRow && cCol)
			 {
			 matrixPrototype[iCell].nullCell = false;

			 matrixPrototype[iCell].reactParams.push_back(
			 r.reactionParametersGet());

			 std::size_t numberOfIndivComps =
			 r.reagentSubstNamesGet().size();
			 if (std::find(r.reagentSubstNamesGet().cbegin(),
			 r.reagentSubstNamesGet().cend(), "M0")
			 != r.reagentSubstNamesGet().cend())
			 numberOfIndivComps--;

			 matrixPrototype[iCell].reactWeight.push_back(
			 1. / scalar(numberOfIndivComps));

			 std::vector<std::pair<std::size_t, int>> compIndexes;
			 for (std::size_t jComp = 0;
			 jComp < r.reagentSubstNamesGet().size();
			 ++jComp)
			 {
			 const auto & rComp = r.reagentSubstNamesGet()[jComp];
			 if (rComp == "M0")
			 compIndexes.emplace_back(
			 std::make_pair(std::size_t(0),
			 r.reagentCoeffsGet()[jComp]));
			 else if (rComp == compNameC)
			 {
			 const auto compC = r.reagentCoeffsGet()[jComp];
			 if (compC != 1)
			 compIndexes.emplace_back(
			 std::make_pair(
			 reactingComponentsMatching.at(
			 rComp) + 1,
			 compC - 1));
			 }
			 else
			 compIndexes.emplace_back(
			 std::make_pair(
			 reactingComponentsMatching.at(
			 rComp) + 1,
			 r.reagentCoeffsGet()[jComp]));
			 }
			 matrixPrototype[iCell].comp.push_back(compIndexes);

			 matrixPrototype[iCell].molMass.second.push_back(
			 therm.Mv()[nR]);

			 for (std::size_t jComp = 0;
			 jComp < r.productSubstNamesGet().size();
			 ++jComp)
			 if (r.productSubstNamesGet()[jComp] == compNameR)
			 {
			 matrixPrototype[iCell].coeffs.push_back(
			 -r.productCoeffsGet()[jComp]);
			 break;
			 }
			 }
			 }
			 }*/
		}

		for (std::size_t iCell = 0; iCell < matrixPrototype.size(); ++iCell)
		{
			std::cout << "Matrix cell: " << iCell << ". Components: "
					<< std::get<0>(matrixPrototype[iCell].indexes) << ' '
					<< std::get<1>(matrixPrototype[iCell].indexes) << '.'
					<< std::endl;
			std::cout << std::boolalpha;
			std::cout << "Null cell: " << matrixPrototype[iCell].nullCell << '.'
					<< std::endl;
			std::cout << "Diagonal cell: "
					<< matrixPrototype[iCell].diagonalCell << '.' << std::endl;
			std::cout << std::noboolalpha;
			std::cout << "Molar mass derivative: "
					<< matrixPrototype[iCell].molMass.first
					<< ", molar mass for reactions: ";
			for (const auto & rM : matrixPrototype[iCell].molMass.second)
				std::cout << rM << ' ';
			std::cout << '.' << std::endl;
			std::cout << "Number of reactions in cell: "
					<< matrixPrototype[iCell].reactParams.size() << std::endl;

			for (std::size_t jReac = 0;
					jReac < matrixPrototype[iCell].reactParams.size(); ++jReac)
			{
				std::cout << "Reaction " << jReac + 1 << '.' << std::endl;
				std::cout << '\t' << "Reaction parameters: " << "A = "
						<< std::get<0>(
								matrixPrototype[iCell].reactParams[jReac])
						<< ", " << "n = "
						<< std::get<1>(
								matrixPrototype[iCell].reactParams[jReac])
						<< ", " << "E = "
						<< std::get<2>(
								matrixPrototype[iCell].reactParams[jReac])
						<< '.' << std::endl;
				std::cout << '\t' << "Reaction weight: "
						<< matrixPrototype[iCell].reactWeight[jReac] << '.'
						<< std::endl;
				std::cout << '\t' << "Reagents and there's exponents: "
						<< std::endl;
				for (std::size_t kComp = 0;
						kComp < matrixPrototype[iCell].comp[jReac].size();
						++kComp)
					std::cout << '\t' << '\t' << "Substance "
							<< matrixPrototype[iCell].comp[jReac][kComp].first
							<< ", it's exponent: "
							<< matrixPrototype[iCell].comp[jReac][kComp].second
							<< '.' << std::endl;

				std::cout << '\t' << "Coefficient of change: "
						<< matrixPrototype[iCell].coeffs[jReac] << '.'
						<< std::endl;
			}
		}
	}
}
