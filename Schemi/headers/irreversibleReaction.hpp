/*
 * reversibleReaction.hpp
 *
 *  Created on: 2026/06/25
 *      Author: Maxim Boldyrev
 */

#ifndef IRREVERSIBLEREACTION_HPP_
#define IRREVERSIBLEREACTION_HPP_

#include <array>
#include <string>
#include <vector>

#include "abstractMixtureThermodynamics.hpp"
#include "scalar.hpp"

namespace schemi
{
namespace chemicalKinetics
{

class irreversibleReaction
{
	std::vector<int> reagentCoeffs {}, productCoeffs {};
	std::vector<std::string> reagentSubstNames {}, productSubstNames {};
	std::array<scalar, 3> reactionParameters { 0, 0, 0 };

public:
	irreversibleReaction(const std::string & reactionDescription);

	const std::vector<int>& reagentCoeffsGet() const noexcept
	{
		return reagentCoeffs;
	}
	const std::vector<int>& productCoeffsGet() const noexcept
	{
		return productCoeffs;
	}
	const std::vector<std::string>& reagentSubstNamesGet() const noexcept
	{
		return reagentSubstNames;
	}
	const std::vector<std::string>& productSubstNamesGet() const noexcept
	{
		return productSubstNames;
	}
	const std::array<scalar, 3>& reactionParametersGet() const noexcept
	{
		return reactionParameters;
	}
};

}
}

#endif /* IRREVERSIBLEREACTION_HPP_ */
