/*
 * typeOfTVDLimiterEnum.hpp
 *
 *  Created on: 2019/12/06
 *      Author: Maxim Boldyrev
 *
 *      Types of TVD-limiter.
 */

#ifndef TYPEOFTVDLIMITERENUM_HPP_
#define TYPEOFTVDLIMITERENUM_HPP_

namespace schemi
{
enum class typeOfTVDLimiterEnum
{
	zero,
	minmod,
	vanLeer,
	linear,
	vanAlbada,
	HQUICK,
	vanLeer2,
	superbee,
	vanAlbada2,
	minmod2
};
}  // namespace schemi

#endif /* TYPEOFTVDLIMITERENUM_HPP_ */
