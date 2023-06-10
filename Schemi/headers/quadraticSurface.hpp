/*
 * quadraticSurface.hpp
 *
 *  Created on: 2019/11/12
 *      Author: Maxim Boldyrev
 *
 *      Structure for quadratic surface.
 */

#ifndef QUADRATICSURFACE_HPP_
#define QUADRATICSURFACE_HPP_

#include "vector.hpp"

namespace schemi
{
struct quadraticSurface
{
	const vector& r00() const noexcept;

	const vector& rX0() const noexcept;

	const vector& r0Y() const noexcept;

	const vector& rXY() const noexcept;

	const vector& rC() const noexcept;

	const scalar& S() const noexcept;

	const vector& N() const noexcept;

	vector& r00_r() noexcept;

	vector& rX0_r() noexcept;

	vector& r0Y_r() noexcept;

	vector& rXY_r() noexcept;

	vector& rC_r() noexcept;

	scalar& S_r() noexcept;

	vector& N_r() noexcept;

	quadraticSurface() noexcept;
private:
	vector rv00, rvX0, rv0Y, rvXY; /*peaks' coordinates*/
	vector rvC; /*coordinates of the center*/
	vector normal; /*normal of surface*/
	scalar Sv = 0.; /*Area of the surface*/
};
}  // namespace schemi

#endif /* QUADRATICSURFACE_HPP_ */
