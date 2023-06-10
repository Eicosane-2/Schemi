/*
 * cubicCell.hpp
 *
 *  Created on: 2019/11/12
 *      Author: Maxim Boldyrev
 *
 *      Class describing cubic cell.
 */

#ifndef CUBICCELL_HPP_
#define CUBICCELL_HPP_

#include "vector.hpp"

namespace schemi
{
struct cubicCell
{
	const vector& r000() const noexcept;

	const vector& rX00() const noexcept;

	const vector& r0Y0() const noexcept;

	const vector& rXY0() const noexcept;

	const vector& r00Z() const noexcept;

	const vector& rX0Z() const noexcept;

	const vector& r0YZ() const noexcept;

	const vector& rXYZ() const noexcept;

	const vector& rC() const noexcept;

	const scalar& V() const noexcept;

	vector& r000_r() noexcept;

	vector& rX00_r() noexcept;

	vector& r0Y0_r() noexcept;

	vector& rXY0_r() noexcept;

	vector& r00Z_r() noexcept;

	vector& rX0Z_r() noexcept;

	vector& r0YZ_r() noexcept;

	vector& rXYZ_r() noexcept;

	vector& rC_r() noexcept;

	scalar& V_r() noexcept;

	cubicCell() noexcept;
private:
	vector rv000, rvX00, rv0Y0, rvXY0, rv00Z, rvX0Z, rv0YZ, rvXYZ; /*peaks' coordinates*/
	vector rvC; /*coordinates of the center*/
	scalar Vv = 0.; /*volume of the cell*/
};
}  // namespace schemi

#endif /* CUBICCELL_HPP_ */
