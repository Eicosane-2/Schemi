/*
 * chemicalReactionsEnum.hpp
 *
 *  Created on: 2023/05/09
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALREACTIONSENUM_HPP_
#define CHEMICALREACTIONSENUM_HPP_

namespace schemi
{
enum class chemicalReactions
{
	noReaction,
	Cl2Dissociation,
	Cl2H2Dissociation,
	H2Cl2Combustion,
	NO2Disproportionation,
	H2O2Combustion
};
}  // namespace schemi

#endif /* CHEMICALREACTIONSENUM_HPP_ */
