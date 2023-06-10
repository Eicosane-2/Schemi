/*
 * enthalpyFlowEnum.hpp
 *
 *  Created on: 2020/07/09
 *      Author: Maxim Boldyrev
 *
 *      Variants of enthalpy flow calculation.
 */

#ifndef ENTHALPYFLOWENUM_HPP_
#define ENTHALPYFLOWENUM_HPP_

namespace schemi
{
enum class enthalpyFlowEnum
{
	explicitSolve, implicitSolve, noSolve
};
}  // namespace schemi

#endif /* ENTHALPYFLOWENUM_HPP_ */
