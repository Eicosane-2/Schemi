/*
 * interfaceStatus.hpp
 *
 *  Created on: 2025/07/11
 *      Author: Maxim Boldyrev
 */

#ifndef INTERFACESTATUSENUM_HPP_
#define INTERFACESTATUSENUM_HPP_

namespace schemi
{
enum class interfaceStatus
{
	notDeveloped, developedResolvable, developedNotResolvable
};

enum class initialisationStatus
{
	notInitialised, initialisedButNotApplied, initialisedAndApplied
};
}  // namespace schemi

#endif /* INTERFACESTATUSENUM_HPP_ */
