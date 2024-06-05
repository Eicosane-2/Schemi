/*
 * chemicalKineticsNoReaction.hpp
 *
 *  Created on: 2023/05/09
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSNOREACTION_HPP_
#define CHEMICALKINETICSNOREACTION_HPP_

#include "abstractChemicalKinetics.hpp"

namespace schemi
{
namespace chemicalKinetics
{
class NoReaction: public abstractChemicalKinetics
{
public:
	NoReaction() noexcept;

	void solveChemicalKinetics(homogeneousPhase<cubicCell>&) const override;
};
}
}  // namespace schemi

#endif /* CHEMICALKINETICSNOREACTION_HPP_ */
