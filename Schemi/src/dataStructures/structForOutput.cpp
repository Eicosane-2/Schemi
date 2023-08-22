/*
 * structForOutput.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "structForOutput.hpp"

schemi::structForOutput::structForOutput(MPIHandler & par,
		const mesh & meshRef_in, const std::size_t numberOfComponents) noexcept :

		x_coord(0), y_coord(0), z_coord(0),

		concentration(numberOfComponents + 1), density(numberOfComponents + 1),

		velocity_x(0), velocity_y(0), velocity_z(0),

		pressure(0), kTurb(0), epsTurb(0),

		aTurb_x(0), aTurb_y(0), aTurb_z(0), bTurb(0),

		momentum_x(0), momentum_y(0), momentum_z(0),

		internalEnergy(0), totalEnergy(0), temperature(0),

		rhokTurb(0), rhoepsTurb(0),

		rhoaTurb_x(0), rhoaTurb_y(0), rhoaTurb_z(0), rhobTurb(0),

		HelmholtzEnergy(0), entropy(0),

		concentrationNonSorted(numberOfComponents), velocity_xNonSorted(0), velocity_yNonSorted(
				0), velocity_zNonSorted(0), pressureNonSorted(0), kTurbNonSorted(
				0), epsTurbNonSorted(0), aTurb_xNonSorted(0), aTurb_yNonSorted(
				0), aTurb_zNonSorted(0), bTurbNonSorted(0),

		tNu(0), sonicSpeed(0),

		parallelismRef(par), meshRef(meshRef_in),

		nodeCellN(
				meshRef.nCells()[0] * meshRef.nCells()[1]
						* meshRef.nCells()[2]),

		localSlice(meshRef.nCells()[0] * meshRef.nCells()[1]),

		globalSlice(
				meshRef.nCells()[0] * parallelismRef.mpi_size
						* meshRef.nCells()[1]),

		k_n(meshRef.nCells()[2]), j_n(meshRef.nCells()[1]), i_n(
				meshRef.nCells()[0] * parallelismRef.mpi_size)
{
}

void schemi::structForOutput::setSizes() noexcept
{
#ifdef MPI_VERSION
	x_coord.resize(parallelismRef.totCellNum());
	y_coord.resize(parallelismRef.totCellNum());
	z_coord.resize(parallelismRef.totCellNum());

	for (auto & c_i : concentration)
		c_i.resize(parallelismRef.totCellNum());

	for (auto & d_i : density)
		d_i.resize(parallelismRef.totCellNum());

	velocity_x.resize(parallelismRef.totCellNum());
	velocity_y.resize(parallelismRef.totCellNum());
	velocity_z.resize(parallelismRef.totCellNum());

	pressure.resize(parallelismRef.totCellNum());

	kTurb.resize(parallelismRef.totCellNum());

	epsTurb.resize(parallelismRef.totCellNum());

	aTurb_x.resize(parallelismRef.totCellNum());
	aTurb_y.resize(parallelismRef.totCellNum());
	aTurb_z.resize(parallelismRef.totCellNum());

	bTurb.resize(parallelismRef.totCellNum());

	momentum_x.resize(parallelismRef.totCellNum());
	momentum_y.resize(parallelismRef.totCellNum());
	momentum_z.resize(parallelismRef.totCellNum());

	internalEnergy.resize(parallelismRef.totCellNum());

	totalEnergy.resize(parallelismRef.totCellNum());

	temperature.resize(parallelismRef.totCellNum());

	rhokTurb.resize(parallelismRef.totCellNum());

	rhoepsTurb.resize(parallelismRef.totCellNum());

	rhoaTurb_x.resize(parallelismRef.totCellNum());
	rhoaTurb_y.resize(parallelismRef.totCellNum());
	rhoaTurb_z.resize(parallelismRef.totCellNum());

	rhobTurb.resize(parallelismRef.totCellNum());

	HelmholtzEnergy.resize(parallelismRef.totCellNum());
	entropy.resize(parallelismRef.totCellNum());

	for (auto & c_i : concentrationNonSorted)
		c_i.resize(parallelismRef.totCellNum());

	velocity_xNonSorted.resize(parallelismRef.totCellNum()), velocity_yNonSorted.resize(
			parallelismRef.totCellNum()), velocity_zNonSorted.resize(
			parallelismRef.totCellNum()), pressureNonSorted.resize(
			parallelismRef.totCellNum()), kTurbNonSorted.resize(
			parallelismRef.totCellNum()), epsTurbNonSorted.resize(
			parallelismRef.totCellNum()), aTurb_xNonSorted.resize(
			parallelismRef.totCellNum()), aTurb_yNonSorted.resize(
			parallelismRef.totCellNum()), aTurb_zNonSorted.resize(
			parallelismRef.totCellNum()), bTurbNonSorted.resize(
			parallelismRef.totCellNum()),

	tNu.resize(parallelismRef.totCellNum());

	sonicSpeed.resize(parallelismRef.totCellNum());
#else
	x_coord.resize(meshRef.cellsSize());
	y_coord.resize(meshRef.cellsSize());
	z_coord.resize(meshRef.cellsSize());

	for (auto & c_i : concentration)
		c_i.resize(meshRef.cellsSize());

	for (auto & d_i : density)
		d_i.resize(meshRef.cellsSize());

	velocity_x.resize(meshRef.cellsSize());
	velocity_y.resize(meshRef.cellsSize());
	velocity_z.resize(meshRef.cellsSize());

	pressure.resize(meshRef.cellsSize());

	kTurb.resize(meshRef.cellsSize());

	epsTurb.resize(meshRef.cellsSize());

	aTurb_x.resize(meshRef.cellsSize());
	aTurb_y.resize(meshRef.cellsSize());
	aTurb_z.resize(meshRef.cellsSize());

	bTurb.resize(meshRef.cellsSize());

	momentum_x.resize(meshRef.cellsSize());
	momentum_y.resize(meshRef.cellsSize());
	momentum_z.resize(meshRef.cellsSize());

	internalEnergy.resize(meshRef.cellsSize());

	totalEnergy.resize(meshRef.cellsSize());

	temperature.resize(meshRef.cellsSize());

	rhokTurb.resize(meshRef.cellsSize());

	rhoepsTurb.resize(meshRef.cellsSize());

	rhoaTurb_x.resize(meshRef.cellsSize());
	rhoaTurb_y.resize(meshRef.cellsSize());
	rhoaTurb_z.resize(meshRef.cellsSize());

	rhobTurb.resize(meshRef.cellsSize());

	HelmholtzEnergy.resize(meshRef.cellsSize());
	entropy.resize(meshRef.cellsSize());

	for (auto & c_i : concentrationNonSorted)
		c_i.resize(meshRef.cellsSize());

	velocity_xNonSorted.resize(meshRef.cellsSize()), velocity_yNonSorted.resize(
			meshRef.cellsSize()), velocity_zNonSorted.resize(
			meshRef.cellsSize()), pressureNonSorted.resize(meshRef.cellsSize()), kTurbNonSorted.resize(
			meshRef.cellsSize()), epsTurbNonSorted.resize(meshRef.cellsSize()), aTurb_xNonSorted.resize(
			meshRef.cellsSize()), aTurb_yNonSorted.resize(meshRef.cellsSize()), aTurb_zNonSorted.resize(
			meshRef.cellsSize()), bTurbNonSorted.resize(meshRef.cellsSize()),

	tNu.resize(meshRef.cellsSize());

	sonicSpeed.resize(meshRef.cellsSize());
#endif
}

void schemi::structForOutput::collectParallelData(
		const bunchOfFields<cubicCell> & cellFields,
		const volumeField<scalar> & tNu_in,
		const volumeField<scalar> & sonicSpeed_in) noexcept
{
	volumeField<scalar> z_node(meshRef, 0);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		z_node.ref_r()[i] = meshRef.cells()[i].rC().v()[2];

	parallelismRef.gatherField(z_node, z_coord);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		z_node.ref_r()[i] = cellFields.velocity.ref()[i].v()[2];

	parallelismRef.gatherField(z_node, velocity_z);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		z_node.ref_r()[i] = cellFields.aTurb.ref()[i].v()[2];

	parallelismRef.gatherField(z_node, aTurb_z);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		z_node.ref_r()[i] = cellFields.momentum.ref()[i].v()[2];

	parallelismRef.gatherField(z_node, momentum_z);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		z_node.ref_r()[i] = cellFields.rhoaTurb.ref()[i].v()[2];

	parallelismRef.gatherField(z_node, rhoaTurb_z);

	volumeField<scalar> y_node(meshRef, 0);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.ref_r()[i] = meshRef.cells()[i].rC().v()[1];

	parallelismRef.gatherField(y_node, y_coord);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.ref_r()[i] = cellFields.velocity.ref()[i].v()[1];

	parallelismRef.gatherField(y_node, velocity_y);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.ref_r()[i] = cellFields.aTurb.ref()[i].v()[1];

	parallelismRef.gatherField(y_node, aTurb_y);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.ref_r()[i] = cellFields.momentum.ref()[i].v()[1];

	parallelismRef.gatherField(y_node, momentum_y);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.ref_r()[i] = cellFields.rhoaTurb.ref()[i].v()[1];

	parallelismRef.gatherField(y_node, rhoaTurb_y);

	volumeField<scalar> x_node(meshRef, 0);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.ref_r()[i] = meshRef.cells()[i].rC().v()[0];

	parallelismRef.gatherField(x_node, x_coord);

	for (std::size_t k = 0; k < cellFields.concentration.v.size(); ++k)
		parallelismRef.gatherField(cellFields.concentration.v[k],
				concentration[k]);

	for (std::size_t k = 0; k < cellFields.density.size(); ++k)
		parallelismRef.gatherField(cellFields.density[k], density[k]);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.ref_r()[i] = cellFields.velocity.ref()[i].v()[0];

	parallelismRef.gatherField(x_node, velocity_x);

	parallelismRef.gatherField(cellFields.pressure, pressure);

	parallelismRef.gatherField(cellFields.kTurb, kTurb);

	parallelismRef.gatherField(cellFields.epsTurb, epsTurb);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.ref_r()[i] = cellFields.aTurb.ref()[i].v()[0];

	parallelismRef.gatherField(x_node, aTurb_x);

	parallelismRef.gatherField(cellFields.bTurb, bTurb);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.ref_r()[i] = cellFields.momentum.ref()[i].v()[0];

	parallelismRef.gatherField(x_node, momentum_x);

	parallelismRef.gatherField(cellFields.internalEnergy, internalEnergy);

	parallelismRef.gatherField(cellFields.totalEnergy, totalEnergy);

	parallelismRef.gatherField(cellFields.temperature, temperature);

	parallelismRef.gatherField(cellFields.rhokTurb, rhokTurb);

	parallelismRef.gatherField(cellFields.rhoepsTurb, rhoepsTurb);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.ref_r()[i] = cellFields.rhoaTurb.ref()[i].v()[0];

	parallelismRef.gatherField(x_node, rhoaTurb_x);

	parallelismRef.gatherField(cellFields.rhobTurb, rhobTurb);

	parallelismRef.gatherField(cellFields.HelmholtzEnergy, HelmholtzEnergy);

	parallelismRef.gatherField(cellFields.entropy, entropy);

	parallelismRef.gatherField(tNu_in, tNu);

	parallelismRef.gatherField(sonicSpeed_in, sonicSpeed);

	if (parallelismRef.isRoot())
	{
		for (std::size_t k = 1; k < concentration.size(); ++k)
			concentrationNonSorted[k - 1] = concentration[k];

		velocity_xNonSorted = velocity_x;
		velocity_yNonSorted = velocity_y;
		velocity_zNonSorted = velocity_z;

		pressureNonSorted = pressure;

		kTurbNonSorted = kTurb;
		epsTurbNonSorted = epsTurb;

		aTurb_xNonSorted = aTurb_x;
		aTurb_yNonSorted = aTurb_y;
		aTurb_zNonSorted = aTurb_z;

		bTurbNonSorted = bTurb;

#ifdef MPI_VERSION
		switch (meshRef.taskDimension())
		{
		case dimensions::task3D:
		case dimensions::task2D:
			rearrange();
			break;
		case dimensions::task1D:
		default:
			break;
		}
#endif
	}
}

void schemi::structForOutput::rearrange() noexcept
{
	const auto x_coord_copy = x_coord;
	const auto y_coord_copy = y_coord;
	const auto z_coord_copy = z_coord;

	const auto concentration_copy = concentration;
	const auto density_copy = density;
	const auto velocity_x_copy = velocity_x;
	const auto velocity_y_copy = velocity_y;
	const auto velocity_z_copy = velocity_z;
	const auto pressure_copy = pressure;
	const auto kTurb_copy = kTurb;
	const auto epsTurb_copy = epsTurb;
	const auto aTurb_x_copy = aTurb_x;
	const auto aTurb_y_copy = aTurb_y;
	const auto aTurb_z_copy = aTurb_z;
	const auto bTurb_copy = bTurb;
	const auto momentum_x_copy = momentum_x;
	const auto momentum_y_copy = momentum_y;
	const auto momentum_z_copy = momentum_z;
	const auto internalEnergy_copy = internalEnergy;
	const auto totalEnergy_copy = totalEnergy;
	const auto temperature_copy = temperature;
	const auto rhokTurb_copy = rhokTurb;
	const auto rhoepsTurb_copy = rhoepsTurb;
	const auto rhoaTurb_x_copy = rhoaTurb_x;
	const auto rhoaTurb_y_copy = rhoaTurb_y;
	const auto rhoaTurb_z_copy = rhoaTurb_z;
	const auto rhobTurb_copy = rhobTurb;
	const auto HelmholtzEnergy_copy = HelmholtzEnergy;
	const auto entropy_copy = entropy;

	const auto tNu_copy = tNu;

	const auto sonicSpeed_copy = sonicSpeed;

	for (std::size_t k = 0; k < k_n; ++k)
		for (std::size_t j = 0; j < j_n; ++j)
			for (std::size_t i = 0; i < i_n; ++i)
			{
				const std::size_t curRank = i / meshRef.nCells()[0];
				const std::size_t rank_i = i % meshRef.nCells()[0];
				const std::size_t prev_cells = curRank * nodeCellN;

				const std::size_t copyFromIndex = prev_cells + k * localSlice
						+ j * meshRef.nCells()[0] + rank_i;
				const std::size_t copyToIndex = i + j * i_n + k * globalSlice;

				rearrange_cell(x_coord, x_coord_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(y_coord, y_coord_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(z_coord, z_coord_copy, copyToIndex,
						copyFromIndex);

				for (std::size_t l = 0; l < concentration.size(); ++l)
				{
					rearrange_cell(concentration[l], concentration_copy[l],
							copyToIndex, copyFromIndex);
					rearrange_cell(density[l], density_copy[l], copyToIndex,
							copyFromIndex);
				}
				rearrange_cell(velocity_x, velocity_x_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(velocity_y, velocity_y_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(velocity_z, velocity_z_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(pressure, pressure_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(kTurb, kTurb_copy, copyToIndex, copyFromIndex);
				rearrange_cell(epsTurb, epsTurb_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(aTurb_x, aTurb_x_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(aTurb_y, aTurb_y_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(aTurb_z, aTurb_z_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(bTurb, bTurb_copy, copyToIndex, copyFromIndex);
				rearrange_cell(momentum_x, momentum_x_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(momentum_y, momentum_y_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(momentum_z, momentum_z_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(internalEnergy, internalEnergy_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(totalEnergy, totalEnergy_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(temperature, temperature_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(rhokTurb, rhokTurb_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(rhoepsTurb, rhoepsTurb_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(rhoaTurb_x, rhoaTurb_x_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(rhoaTurb_y, rhoaTurb_y_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(rhoaTurb_z, rhoaTurb_z_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(rhobTurb, rhobTurb_copy, copyToIndex,
						copyFromIndex);
				rearrange_cell(HelmholtzEnergy, HelmholtzEnergy_copy,
						copyToIndex, copyFromIndex);
				rearrange_cell(entropy, entropy_copy, copyToIndex,
						copyFromIndex);

				rearrange_cell(tNu, tNu_copy, copyToIndex, copyFromIndex);

				rearrange_cell(sonicSpeed, sonicSpeed_copy, copyToIndex,
						copyFromIndex);
			}
}

void schemi::structForOutput::rearrange_cell(std::valarray<scalar> & to,
		const std::valarray<scalar> & from, const std::size_t copyToIndex,
		const std::size_t copyFromIndex) const noexcept
{
	to[copyToIndex] = from[copyFromIndex];
}
