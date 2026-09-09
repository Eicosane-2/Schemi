/*
 * reversibleReaction.cpp
 *
 *  Created on: 2026/06/25
 *      Author: Maxim Boldyrev
 */

#include "irreversibleReaction.hpp"

#include <algorithm>
#include <sstream>

#include "split.hpp"

schemi::chemicalKinetics::irreversibleReaction::irreversibleReaction(
		const std::string & reactionDescription)
{
	const auto [reaction, parameters] = split(reactionDescription, '|');

	const auto [forwReaction, products] = split(reaction, '=');

	std::istringstream paramStream(parameters);

	paramStream >> std::get<0>(reactionParameters)
			>> std::get<1>(reactionParameters)
			>> std::get<2>(reactionParameters);

	const std::size_t numberOfCompReag = std::count(forwReaction.cbegin(),
			forwReaction.cend(), '+') + 1;
	const std::size_t numberOfCompProd = std::count(products.cbegin(),
			products.cend(), '+') + 1;

	reagentCoeffs.resize(numberOfCompReag), reagentSubstNames.resize(
			numberOfCompReag);
	productCoeffs.resize(numberOfCompProd), productSubstNames.resize(
			numberOfCompProd);

	std::vector<std::string::const_iterator> plusPositionsReag,
			plusPositionsProd;

	if (numberOfCompReag == 1)
	{
		const auto & coeffComp = forwReaction;

		const auto [coeff, comp] = split(coeffComp, '*');

		reagentCoeffs[0] = std::stoi(coeff);
		reagentSubstNames[0] = comp;
	}
	else
	{
		for (auto it = forwReaction.cbegin(); it != forwReaction.cend(); ++it)
			if (*it == '+')
				plusPositionsReag.push_back(it);

		for (std::size_t k = 0; k < numberOfCompReag; k++)
		{
			if (k == 0)
			{
				const auto coeffComp = std::string(forwReaction.cbegin(),
						plusPositionsReag[k]);

				const auto [coeff, comp] = split(coeffComp, '*');

				reagentCoeffs[k] = std::stoi(coeff);
				reagentSubstNames[k] = comp;
			}
			else if (k == numberOfCompReag - 1)
			{
				auto beg = plusPositionsReag[plusPositionsReag.size() - 1];
				beg++;
				const auto coeffComp = std::string(beg, forwReaction.cend());

				const auto [coeff, comp] = split(coeffComp, '*');

				reagentCoeffs[k] = std::stoi(coeff);
				reagentSubstNames[k] = comp;
			}
			else
			{
				auto beg = plusPositionsReag[k - 1];
				beg++;
				const auto coeffComp = std::string(beg, plusPositionsReag[k]);

				const auto [coeff, comp] = split(coeffComp, '*');

				reagentCoeffs[k] = std::stoi(coeff);
				reagentSubstNames[k] = comp;
			}
		}
	}

	if (numberOfCompProd == 1)
	{
		const auto & coeffComp = products;

		const auto [coeff, comp] = split(coeffComp, '*');

		productCoeffs[0] = std::stoi(coeff);
		productSubstNames[0] = comp;
	}
	else
	{
		for (auto it = products.cbegin(); it != products.cend(); ++it)
			if (*it == '+')
				plusPositionsProd.push_back(it);

		for (std::size_t k = 0; k < numberOfCompProd; k++)
		{
			if (k == 0)
			{
				const auto coeffComp = std::string(products.cbegin(),
						plusPositionsProd[k]);

				const auto [coeff, comp] = split(coeffComp, '*');

				productCoeffs[k] = std::stoi(coeff);
				productSubstNames[k] = comp;
			}
			else if (k == numberOfCompProd - 1)
			{
				auto beg = plusPositionsProd[plusPositionsProd.size() - 1];
				beg++;
				const auto coeffComp = std::string(beg, products.cend());

				const auto [coeff, comp] = split(coeffComp, '*');

				productCoeffs[k] = std::stoi(coeff);
				productSubstNames[k] = comp;
			}
			else
			{
				auto beg = plusPositionsProd[k - 1];
				beg++;
				const auto coeffComp = std::string(beg, plusPositionsProd[k]);

				const auto [coeff, comp] = split(coeffComp, '*');

				productCoeffs[k] = std::stoi(coeff);
				productSubstNames[k] = comp;
			}
		}
	}
}
