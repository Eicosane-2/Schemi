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

		parallelismRef(par), meshRef(meshRef_in)
{
}

void schemi::structForOutput::setSizes() noexcept
{
#ifdef MPI_VERSION
	if (parallelismRef.isRoot())
	{
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
	}
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
		z_node.r()[i] = std::get<2>(meshRef.cells()[i].rC()());

	parallelismRef.gatherField(z_node, z_coord);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		z_node.r()[i] = std::get<2>(cellFields.velocity()[i]());

	parallelismRef.gatherField(z_node, velocity_z);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		z_node.r()[i] = std::get<2>(cellFields.aTurb()[i]());

	parallelismRef.gatherField(z_node, aTurb_z);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		z_node.r()[i] = std::get<2>(cellFields.momentum()[i]());

	parallelismRef.gatherField(z_node, momentum_z);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		z_node.r()[i] = std::get<2>(cellFields.rhoaTurb()[i]());

	parallelismRef.gatherField(z_node, rhoaTurb_z);

	volumeField<scalar> y_node(meshRef, 0);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.r()[i] = std::get<1>(meshRef.cells()[i].rC()());

	parallelismRef.gatherField(y_node, y_coord);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.r()[i] = std::get<1>(cellFields.velocity()[i]());

	parallelismRef.gatherField(y_node, velocity_y);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.r()[i] = std::get<1>(cellFields.aTurb()[i]());

	parallelismRef.gatherField(y_node, aTurb_y);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.r()[i] = std::get<1>(cellFields.momentum()[i]());

	parallelismRef.gatherField(y_node, momentum_y);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		y_node.r()[i] = std::get<1>(cellFields.rhoaTurb()[i]());

	parallelismRef.gatherField(y_node, rhoaTurb_y);

	volumeField<scalar> x_node(meshRef, 0);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.r()[i] = std::get<0>(meshRef.cells()[i].rC()());

	parallelismRef.gatherField(x_node, x_coord);

	for (std::size_t k = 0; k < cellFields.concentration.v.size(); ++k)
		parallelismRef.gatherField(cellFields.concentration.v[k],
				concentration[k]);

	for (std::size_t k = 0; k < cellFields.density.size(); ++k)
		parallelismRef.gatherField(cellFields.density[k], density[k]);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.r()[i] = std::get<0>(cellFields.velocity()[i]());

	parallelismRef.gatherField(x_node, velocity_x);

	parallelismRef.gatherField(cellFields.pressure, pressure);

	parallelismRef.gatherField(cellFields.kTurb, kTurb);

	parallelismRef.gatherField(cellFields.epsTurb, epsTurb);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.r()[i] = std::get<0>(cellFields.aTurb()[i]());

	parallelismRef.gatherField(x_node, aTurb_x);

	parallelismRef.gatherField(cellFields.bTurb, bTurb);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.r()[i] = std::get<0>(cellFields.momentum()[i]());

	parallelismRef.gatherField(x_node, momentum_x);

	parallelismRef.gatherField(cellFields.internalEnergy, internalEnergy);

	parallelismRef.gatherField(cellFields.totalEnergy, totalEnergy);

	parallelismRef.gatherField(cellFields.temperature, temperature);

	parallelismRef.gatherField(cellFields.rhokTurb, rhokTurb);

	parallelismRef.gatherField(cellFields.rhoepsTurb, rhoepsTurb);

	for (std::size_t i = 0; i < meshRef.cellsSize(); ++i)
		x_node.r()[i] = std::get<0>(cellFields.rhoaTurb()[i]());

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
			rearrange3D();
			break;
		case dimensions::task2D:
			rearrange2D();
			break;
		case dimensions::task1D:
		default:
			break;
		}
#endif
	}
}

void schemi::structForOutput::rearrange2D() noexcept
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

	const std::size_t cellsPerNode = meshRef.cellsSize();
	const std::size_t xMPI_length = std::get<0>(
			parallelismRef.rankDirectionMax());
	const std::size_t cellsPerXNode = std::get<0>(meshRef.nCells());
	const std::size_t cellsPerYNode = std::get<1>(meshRef.nCells());
	const std::size_t cellPerXGlobal = cellsPerXNode * xMPI_length;

	for (std::size_t i = 0; i < x_coord.size(); ++i)
	{
		const std::size_t rank = i / cellsPerNode;
		const std::size_t yRank = rank / xMPI_length;
		const std::size_t xRank = rank % xMPI_length;

		const std::size_t local_i = i % cellsPerNode;
		const std::size_t local_y = local_i / cellsPerXNode;
		const std::size_t local_x = local_i % cellsPerXNode;

		const std::size_t newX = xRank * cellsPerXNode + local_x;
		const std::size_t newY = yRank * cellsPerYNode + local_y;
		const std::size_t newIndex = newY * cellPerXGlobal + newX;

		copy_cell(x_coord, x_coord_copy, newIndex, i);
		copy_cell(y_coord, y_coord_copy, newIndex, i);
		copy_cell(z_coord, z_coord_copy, newIndex, i);

		for (std::size_t l = 0; l < concentration.size(); ++l)
		{
			copy_cell(concentration[l], concentration_copy[l], newIndex, i);
			copy_cell(density[l], density_copy[l], newIndex, i);
		}

		copy_cell(velocity_x, velocity_x_copy, newIndex, i);
		copy_cell(velocity_y, velocity_y_copy, newIndex, i);
		copy_cell(velocity_z, velocity_z_copy, newIndex, i);
		copy_cell(pressure, pressure_copy, newIndex, i);
		copy_cell(kTurb, kTurb_copy, newIndex, i);
		copy_cell(epsTurb, epsTurb_copy, newIndex, i);
		copy_cell(aTurb_x, aTurb_x_copy, newIndex, i);
		copy_cell(aTurb_y, aTurb_y_copy, newIndex, i);
		copy_cell(aTurb_z, aTurb_z_copy, newIndex, i);
		copy_cell(bTurb, bTurb_copy, newIndex, i);
		copy_cell(momentum_x, momentum_x_copy, newIndex, i);
		copy_cell(momentum_y, momentum_y_copy, newIndex, i);
		copy_cell(momentum_z, momentum_z_copy, newIndex, i);
		copy_cell(internalEnergy, internalEnergy_copy, newIndex, i);
		copy_cell(totalEnergy, totalEnergy_copy, newIndex, i);
		copy_cell(temperature, temperature_copy, newIndex, i);
		copy_cell(rhokTurb, rhokTurb_copy, newIndex, i);
		copy_cell(rhoepsTurb, rhoepsTurb_copy, newIndex, i);
		copy_cell(rhoaTurb_x, rhoaTurb_x_copy, newIndex, i);
		copy_cell(rhoaTurb_y, rhoaTurb_y_copy, newIndex, i);
		copy_cell(rhoaTurb_z, rhoaTurb_z_copy, newIndex, i);
		copy_cell(rhobTurb, rhobTurb_copy, newIndex, i);
		copy_cell(HelmholtzEnergy, HelmholtzEnergy_copy, newIndex, i);
		copy_cell(entropy, entropy_copy, newIndex, i);

		copy_cell(tNu, tNu_copy, newIndex, i);

		copy_cell(sonicSpeed, sonicSpeed_copy, newIndex, i);
	}
}

void schemi::structForOutput::rearrange3D() noexcept
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

	const std::size_t cellsPerNode = meshRef.cellsSize();
	const std::size_t xMPI_length = std::get<0>(
			parallelismRef.rankDirectionMax());
	const std::size_t yMPI_length = std::get<1>(
			parallelismRef.rankDirectionMax());
	const std::size_t xyMPI_slice = xMPI_length * yMPI_length;
	const std::size_t cellsPerXNode = std::get<0>(meshRef.nCells());
	const std::size_t cellsPerYNode = std::get<1>(meshRef.nCells());
	const std::size_t cellsPerZNode = std::get<2>(meshRef.nCells());
	const std::size_t xySlice = cellsPerXNode * cellsPerYNode;
	const std::size_t cellPerXGlobal = cellsPerXNode * xMPI_length;
	const std::size_t cellPerYGlobal = cellsPerYNode * yMPI_length;

	for (std::size_t i = 0; i < x_coord.size(); ++i)
	{
		const std::size_t rank = i / cellsPerNode;
		const std::size_t zRank = rank / xyMPI_slice;
		const std::size_t yRank = rank % xyMPI_slice / xMPI_length;
		const std::size_t xRank = rank % xyMPI_slice % xMPI_length;

		const std::size_t local_i = i % cellsPerNode;
		const std::size_t local_z = local_i / xySlice;
		const std::size_t local_y = local_i % xySlice / cellsPerXNode;
		const std::size_t local_x = local_i % xySlice % cellsPerXNode;

		const std::size_t newX = xRank * cellsPerXNode + local_x;
		const std::size_t newY = yRank * cellsPerYNode + local_y;
		const std::size_t newZ = zRank * cellsPerZNode + local_z;
		const std::size_t newIndex = newZ * cellPerXGlobal * cellPerYGlobal
				+ newY * cellPerXGlobal + newX;

		copy_cell(x_coord, x_coord_copy, newIndex, i);
		copy_cell(y_coord, y_coord_copy, newIndex, i);
		copy_cell(z_coord, z_coord_copy, newIndex, i);

		for (std::size_t l = 0; l < concentration.size(); ++l)
		{
			copy_cell(concentration[l], concentration_copy[l], newIndex, i);
			copy_cell(density[l], density_copy[l], newIndex, i);
		}

		copy_cell(velocity_x, velocity_x_copy, newIndex, i);
		copy_cell(velocity_y, velocity_y_copy, newIndex, i);
		copy_cell(velocity_z, velocity_z_copy, newIndex, i);
		copy_cell(pressure, pressure_copy, newIndex, i);
		copy_cell(kTurb, kTurb_copy, newIndex, i);
		copy_cell(epsTurb, epsTurb_copy, newIndex, i);
		copy_cell(aTurb_x, aTurb_x_copy, newIndex, i);
		copy_cell(aTurb_y, aTurb_y_copy, newIndex, i);
		copy_cell(aTurb_z, aTurb_z_copy, newIndex, i);
		copy_cell(bTurb, bTurb_copy, newIndex, i);
		copy_cell(momentum_x, momentum_x_copy, newIndex, i);
		copy_cell(momentum_y, momentum_y_copy, newIndex, i);
		copy_cell(momentum_z, momentum_z_copy, newIndex, i);
		copy_cell(internalEnergy, internalEnergy_copy, newIndex, i);
		copy_cell(totalEnergy, totalEnergy_copy, newIndex, i);
		copy_cell(temperature, temperature_copy, newIndex, i);
		copy_cell(rhokTurb, rhokTurb_copy, newIndex, i);
		copy_cell(rhoepsTurb, rhoepsTurb_copy, newIndex, i);
		copy_cell(rhoaTurb_x, rhoaTurb_x_copy, newIndex, i);
		copy_cell(rhoaTurb_y, rhoaTurb_y_copy, newIndex, i);
		copy_cell(rhoaTurb_z, rhoaTurb_z_copy, newIndex, i);
		copy_cell(rhobTurb, rhobTurb_copy, newIndex, i);
		copy_cell(HelmholtzEnergy, HelmholtzEnergy_copy, newIndex, i);
		copy_cell(entropy, entropy_copy, newIndex, i);

		copy_cell(tNu, tNu_copy, newIndex, i);

		copy_cell(sonicSpeed, sonicSpeed_copy, newIndex, i);
	}
}

void schemi::structForOutput::copy_cell(std::valarray<scalar> & to,
		const std::valarray<scalar> & from, const std::size_t copyToIndex,
		const std::size_t copyFromIndex) const noexcept
{
	to[copyToIndex] = from[copyFromIndex];
}
