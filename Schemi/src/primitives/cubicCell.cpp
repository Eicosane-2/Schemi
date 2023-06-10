/*
 * cubicCell.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "cubicCell.hpp"

const schemi::vector& schemi::cubicCell::r000() const noexcept
{
	return rv000;
}

const schemi::vector& schemi::cubicCell::rX00() const noexcept
{
	return rvX00;
}

const schemi::vector& schemi::cubicCell::r0Y0() const noexcept
{
	return rv0Y0;
}

const schemi::vector& schemi::cubicCell::rXY0() const noexcept
{
	return rvXY0;
}

const schemi::vector& schemi::cubicCell::r00Z() const noexcept
{
	return rv00Z;
}

const schemi::vector& schemi::cubicCell::rX0Z() const noexcept
{
	return rvX0Z;
}

const schemi::vector& schemi::cubicCell::r0YZ() const noexcept
{
	return rv0YZ;
}

const schemi::vector& schemi::cubicCell::rXYZ() const noexcept
{
	return rvXYZ;
}

const schemi::vector& schemi::cubicCell::rC() const noexcept
{
	return rvC;
}

const schemi::scalar& schemi::cubicCell::V() const noexcept
{
	return Vv;
}

schemi::vector& schemi::cubicCell::r000_r() noexcept
{
	return rv000;
}

schemi::vector& schemi::cubicCell::rX00_r() noexcept
{
	return rvX00;
}

schemi::vector& schemi::cubicCell::r0Y0_r() noexcept
{
	return rv0Y0;
}

schemi::vector& schemi::cubicCell::rXY0_r() noexcept
{
	return rvXY0;
}

schemi::vector& schemi::cubicCell::r00Z_r() noexcept
{
	return rv00Z;
}

schemi::vector& schemi::cubicCell::rX0Z_r() noexcept
{
	return rvX0Z;
}

schemi::vector& schemi::cubicCell::r0YZ_r() noexcept
{
	return rv0YZ;
}

schemi::vector& schemi::cubicCell::rXYZ_r() noexcept
{
	return rvXYZ;
}

schemi::vector& schemi::cubicCell::rC_r() noexcept
{
	return rvC;
}

schemi::scalar& schemi::cubicCell::V_r() noexcept
{
	return Vv;
}

schemi::cubicCell::cubicCell() noexcept :
		rv000(0), rvX00(0), rv0Y0(0), rvXY0(0), rv00Z(0), rvX0Z(0), rv0YZ(0), rvXYZ(
				0), rvC(0)
{
}
