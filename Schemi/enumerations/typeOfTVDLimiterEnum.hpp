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
enum class typeOfTVDLimiter
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
	minmod2,
	Sweby
};
}  // namespace schemi

#endif /* TYPEOFTVDLIMITERENUM_HPP_ */
