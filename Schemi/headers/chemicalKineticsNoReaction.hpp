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
class chemicalKineticsNoReaction: public abstractChemicalKinetics
{
public:
	chemicalKineticsNoReaction() noexcept;

	void solveChemicalKinetics(homogeneousPhase<cubicCell>&) const override;
};
}  // namespace schemi

#endif /* CHEMICALKINETICSNOREACTION_HPP_ */
