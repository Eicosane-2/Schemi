/*
 * gasModelEnum.hpp
 *
 *  Created on: 2019/12/01
 *      Author: Maxim Boldyrev
 *
 *      Types of equation of state.
 */

#ifndef GASMODELENUM_HPP_
#define GASMODELENUM_HPP_

namespace schemi
{
enum class gasModel
{
	ideal, vanDerWaals, RedlichKwong, stiffened
};
}  // namespace schemi

#endif /* GASMODELENUM_HPP_ */
