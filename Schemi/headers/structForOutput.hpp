/*
 * structForOutput.hpp
 *
 *  Created on: 2020/09/22
 *      Author: Maxim Boldyrev
 *
 *      Structure where data for output collected.
 */

#ifndef STRUCTFOROUTPUT_HPP_
#define STRUCTFOROUTPUT_HPP_

#include <valarray>
#include <vector>

#include "mesh.hpp"
#include "MPIHandler.hpp"
#include "scalar.hpp"
#include "bunchOfFields.hpp"
#include "volumeField.hpp"

namespace schemi
{
struct structForOutput
{
	structForOutput(const structForOutput&) = delete;
	auto& operator=(const structForOutput&) = delete;

	std::valarray<scalar> x_coord;
	std::valarray<scalar> y_coord;
	std::valarray<scalar> z_coord;

	std::vector<std::valarray<scalar>> concentration;
	std::vector<std::valarray<scalar>> density;
	std::valarray<scalar> velocity_x;
	std::valarray<scalar> velocity_y;
	std::valarray<scalar> velocity_z;
	std::valarray<scalar> pressure;
	std::valarray<scalar> kTurb;
	std::valarray<scalar> epsTurb;
	std::valarray<scalar> aTurb_x;
	std::valarray<scalar> aTurb_y;
	std::valarray<scalar> aTurb_z;
	std::valarray<scalar> bTurb;
	std::valarray<scalar> momentum_x;
	std::valarray<scalar> momentum_y;
	std::valarray<scalar> momentum_z;
	std::valarray<scalar> internalEnergy;
	std::valarray<scalar> totalEnergy;
	std::valarray<scalar> temperature;
	std::valarray<scalar> rhokTurb;
	std::valarray<scalar> rhoepsTurb;
	std::valarray<scalar> rhoaTurb_x;
	std::valarray<scalar> rhoaTurb_y;
	std::valarray<scalar> rhoaTurb_z;
	std::valarray<scalar> rhobTurb;
	std::valarray<scalar> HelmholtzEnergy;
	std::valarray<scalar> entropy;

	std::vector<std::valarray<scalar>> concentrationNonSorted;
	std::valarray<scalar> velocity_xNonSorted;
	std::valarray<scalar> velocity_yNonSorted;
	std::valarray<scalar> velocity_zNonSorted;
	std::valarray<scalar> pressureNonSorted;
	std::valarray<scalar> kTurbNonSorted;
	std::valarray<scalar> epsTurbNonSorted;
	std::valarray<scalar> aTurb_xNonSorted;
	std::valarray<scalar> aTurb_yNonSorted;
	std::valarray<scalar> aTurb_zNonSorted;
	std::valarray<scalar> bTurbNonSorted;

	std::valarray<scalar> tNu;

	std::valarray<scalar> sonicSpeed;

	structForOutput(MPIHandler & par, const mesh & meshRef_in,
			const std::size_t numberOfComponents) noexcept;

	void setSizes() noexcept;

	void collectParallelData(const bunchOfFields<cubicCell> & cellFields,
			const volumeField<scalar> & tNu_in,
			const volumeField<scalar> & sonicSpeed_in) noexcept;
private:
	MPIHandler & parallelismRef;
	const mesh & meshRef;

	const std::size_t nodeCellN, localSlice, globalSlice, k_n, j_n, i_n;

	void rearrange() noexcept;

	void rearrange_cell(std::valarray<scalar> & to,
			const std::valarray<scalar> & from, const std::size_t copyToIndex,
			const std::size_t copyFromIndex) const noexcept;
};
}  // namespace schemi

#endif /* STRUCTFOROUTPUT_HPP_ */
