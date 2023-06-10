/*
 * pressureStarClass.hpp
 *
 *  Created on: 2020/03/31
 *      Author: Maxim Boldyrev
 *
 *      Class for pressure calculation in star zone through equation of state.
 */

#ifndef PRESSURESTARCLASS_HPP_
#define PRESSURESTARCLASS_HPP_

#include "abstractMixtureThermodynamics.hpp"
#include "scalar.hpp"
#include "vector.hpp"

namespace schemi
{
class pressureStarClass
{
protected:
	scalar pressureStar(const abstractMixtureThermodynamics & mix,
			const std::valarray<scalar> & densityState,
			const vector & momentumState, const scalar totalEnergyState,
			const scalar rhokState) const noexcept;
public:
	virtual ~pressureStarClass() noexcept =0;
};
}  // namespace schemi

#endif /* PRESSURESTARCLASS_HPP_ */
