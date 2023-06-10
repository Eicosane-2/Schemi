/*
 * flowSolverEnum.hpp
 *
 *  Created on: 2020/03/31
 *      Author: Maxim Boldyrev
 *
 *      Types of Riemann solver.
 */

#ifndef FLOWSOLVERENUM_HPP_
#define FLOWSOLVERENUM_HPP_

namespace schemi
{
enum class flowSolverEnum
{
	HLL, HLLCF, HLLC, HLLCLM, KT, HLLC2p, Richtmyer
};
}  // namespace schemi

#endif /* FLOWSOLVERENUM_HPP_ */
