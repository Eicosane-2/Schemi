/*
 * mesh.hpp
 *
 *  Created on: 2019/11/12
 *      Author: Maxim Boldyrev
 *
 *      Class for mesh singleton.
 */

#ifndef MESH_HPP_
#define MESH_HPP_

#include <vector>

#include "boundaryConditionTypesEnum.hpp"
#include "dimensionsEnum.hpp"
#include "cubicCell.hpp"
#include "quadraticSurface.hpp"

namespace schemi
{
class mesh;

struct meshDestroyer
{
	meshDestroyer() noexcept;

	~meshDestroyer() noexcept;

	meshDestroyer(const meshDestroyer&) = delete;
	auto& operator=(const meshDestroyer&) = delete;

	void setPointer(mesh * in_p) noexcept;

private:
	mesh * p;
};

class mesh
{
	friend struct meshDestroyer;

	std::vector<cubicCell> cellsA;
	std::vector<quadraticSurface> surfacesA;

	std::vector<std::vector<std::size_t>> surfacesOfCellsA;
	std::vector<std::vector<std::size_t>> neighboursOfCellsA;

	std::vector<std::size_t> surfaceOwnerA;
	std::vector<std::size_t> surfaceNeighbourA;

	std::vector<boundaryConditionType> surfaceBoundaryCondition;

	std::vector<scalar> surfOwnWA;
	std::vector<scalar> surfNeiWA;

	std::vector<std::pair<scalar, std::vector<scalar>>> cellSurfaceDistancesA;
	std::vector<std::pair<scalar, std::vector<scalar>>> cellSurfaceDistancesRA;

	std::size_t tailSurfacesNumber = 0;
	std::size_t innerSurfacesNumber = 0;
	std::size_t pointSurfacesNumber = 0;
	std::size_t bottomSurfacesNumber = 0;
	std::size_t rightSurfacesNumber = 0;
	std::size_t leftSurfacesNumber = 0;
	std::size_t topSurfacesNumber = 0;

	std::size_t nonexistentCell = static_cast<std::size_t>(-1);

	std::size_t cellsNumber = 0;
	std::size_t surfacesNumber = 0;

	bool initialised = false;

	dimensions taskDim = dimensions::task1D;

	std::array<std::size_t, 3> n_cells { 0, 0, 0 };

	vector deltaParallelepiped = vector(0), zeroPoint = vector(0);

	scalar timestep_val, timestepSource_val;

	void calculateNormales() noexcept;

	void calculateWeights() noexcept;

	void calculateCellSurfaceDistances() noexcept;

	mesh() noexcept;

	mesh(const mesh&) = delete;
	auto& operator=(const mesh&) = delete;

	~mesh() noexcept;

	static mesh * pInstance;
	static meshDestroyer destroyer;
public:
	static mesh* instance();

	bool is_initialised() const noexcept;

	const std::vector<cubicCell>& cells() const noexcept;

	const std::vector<quadraticSurface>& surfaces() const noexcept;

	const std::vector<std::size_t>& surfaceOwner() const noexcept;

	const std::vector<std::size_t>& surfaceNeighbour() const noexcept;

	const std::vector<std::vector<std::size_t>>& surfacesOfCells() const noexcept;

	const std::vector<std::vector<std::size_t>>& neighboursOfCells() const noexcept;

	const std::vector<boundaryConditionType>& bndType() const noexcept;

	const std::vector<scalar>& surfOwnW() const noexcept;
	const std::vector<scalar>& surfNeiW() const noexcept;

	const std::vector<std::pair<scalar, std::vector<scalar>>>& cellSurfaceDistances() const noexcept;
	const std::vector<std::pair<scalar, std::vector<scalar>>>& cellSurfaceDistancesR() const noexcept;

	std::size_t tailNumber() const noexcept;

	std::size_t innerNumber() const noexcept;

	std::size_t pointNumber() const noexcept;

	std::size_t bottomNumber() const noexcept;

	std::size_t rightNumber() const noexcept;

	std::size_t leftNumber() const noexcept;

	std::size_t topNumber() const noexcept;

	std::size_t nonexistCell() const noexcept;

	std::size_t cellsSize() const noexcept;

	std::size_t surfacesSize() const noexcept;

	dimensions taskDimension() const noexcept;

	const std::array<std::size_t, 3>& nCells() const noexcept;

	const vector& delta() const noexcept;

	const vector& zero() const noexcept;

	scalar timestep() const noexcept;

	void setTimestep(const scalar dt) noexcept;

	scalar timestepSource() const noexcept;

	scalar& timestepSourceRef() noexcept;

	void oneDParallelepiped(
			const std::pair<vector, vector> & vectorOfParallelepiped,
			const std::size_t N_x,
			const std::vector<boundaryConditionType> & commonConditions);

	void twoDParallelepiped(
			const std::pair<vector, vector> & vectorOfParallelepiped,
			const std::size_t N_x, const std::size_t N_y,
			const std::vector<boundaryConditionType> & commonConditions);

	void threeDParallelepiped(
			const std::pair<vector, vector> & vectorOfParallelepiped,
			const std::size_t N_x, const std::size_t N_y, const std::size_t N_z,
			const std::vector<boundaryConditionType> & commonConditions);

	std::size_t findSeparatingSurface(std::size_t cell1,
			std::size_t cell2) const;
};
}  // namespace schemi

#endif /* MESH_HPP_ */
