/*
 * quadraticSurface.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "quadraticSurface.hpp"

const schemi::vector& schemi::quadraticSurface::r00() const noexcept
{
	return rv00;
}

const schemi::vector& schemi::quadraticSurface::rX0() const noexcept
{
	return rvX0;
}

const schemi::vector& schemi::quadraticSurface::r0Y() const noexcept
{
	return rv0Y;
}

const schemi::vector& schemi::quadraticSurface::rXY() const noexcept
{
	return rvXY;
}

const schemi::vector& schemi::quadraticSurface::rC() const noexcept
{
	return rvC;
}

const schemi::scalar& schemi::quadraticSurface::S() const noexcept
{
	return Sv;
}

const schemi::vector& schemi::quadraticSurface::N() const noexcept
{
	return normal;
}

schemi::vector& schemi::quadraticSurface::r00_r() noexcept
{
	return rv00;
}

schemi::vector& schemi::quadraticSurface::rX0_r() noexcept
{
	return rvX0;
}

schemi::vector& schemi::quadraticSurface::r0Y_r() noexcept
{
	return rv0Y;
}

schemi::vector& schemi::quadraticSurface::rXY_r() noexcept
{
	return rvXY;
}

schemi::vector& schemi::quadraticSurface::rC_r() noexcept
{
	return rvC;
}

schemi::scalar& schemi::quadraticSurface::S_r() noexcept
{
	return Sv;
}

schemi::vector& schemi::quadraticSurface::N_r() noexcept
{
	return normal;
}

schemi::quadraticSurface::quadraticSurface() noexcept :
		rv00(0), rvX0(0), rv0Y(0), rvXY(0), rvC(0), normal(0)
{
}

