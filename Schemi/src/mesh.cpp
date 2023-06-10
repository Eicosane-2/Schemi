/*
 * mesh.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "mesh.hpp"

#include <iostream>

#include "exception.hpp"
#include "globalConstants.hpp"
#include "vector.hpp"
#include "vectorProduct.hpp"
#include "vectorVectorDotProduct.hpp"

schemi::meshDestroyer::meshDestroyer() noexcept :
		p(nullptr)
{
}

void schemi::meshDestroyer::setPointer(mesh * in_p) noexcept
{
	p = in_p;
}

schemi::meshDestroyer::~meshDestroyer() noexcept
{
	delete p;
}

schemi::meshDestroyer schemi::mesh::destroyer;

void schemi::mesh::calculateNormales() noexcept
{
	for (std::size_t i = 0; i < surfacesA.size(); ++i)
	{
		const vector V1(surfacesA[i].rX0() - surfacesA[i].r00());
		const vector V2(surfacesA[i].r0Y() - surfacesA[i].r00());

		vector Normal = vectorProduct(V1, V2);
		Normal /= Normal.mag();

		const std::size_t cellIndex = surfaceOwnerA[i];
		const vector cellSurfaceVector = surfacesA[i].rC()
				- cellsA[cellIndex].rC();

		if ((cellSurfaceVector & Normal) < 0.)
			Normal *= (-1.);

		surfacesA[i].N_r() = Normal;
	}
}

schemi::mesh::mesh() noexcept :
		cellsA(0), surfacesA(0), surfacesOfCellsA(0), neighboursOfCellsA(0), surfaceOwnerA(
				0), surfaceNeighbourA(0), surfaceBoundaryCondition(0), timestep_val(
				veryBig), timestepSource_val(veryBig)
{
}

schemi::mesh::~mesh() noexcept
{
	std::cout << "Mesh destroyed." << std::endl;
}

schemi::mesh* schemi::mesh::instance()
{
	if (pInstance == nullptr)
	{
		pInstance = new mesh();

		destroyer.setPointer(pInstance);

		return pInstance;
	}
	else
		throw exception("Mesh object should be created only once.",
				errorsEnum::meshGenerationError);
}

bool schemi::mesh::is_initialized() const noexcept
{
	return initialized;
}

const std::vector<schemi::cubicCell>& schemi::mesh::cells() const noexcept
{
	return cellsA;
}

const std::vector<schemi::quadraticSurface>& schemi::mesh::surfaces() const noexcept
{
	return surfacesA;
}

const std::vector<std::size_t>& schemi::mesh::surfaceOwner() const noexcept
{
	return surfaceOwnerA;
}

const std::vector<std::size_t>& schemi::mesh::surfaceNeighbour() const noexcept
{
	return surfaceNeighbourA;
}

const std::vector<std::vector<std::size_t>>& schemi::mesh::surfacesOfCells() const noexcept
{
	return surfacesOfCellsA;
}

const std::vector<std::vector<std::size_t>>& schemi::mesh::neighboursOfCells() const noexcept
{
	return neighboursOfCellsA;
}

const std::vector<schemi::boundaryConditionType>& schemi::mesh::bndType() const noexcept
{
	return surfaceBoundaryCondition;
}

std::size_t schemi::mesh::tailNumber() const noexcept
{
	return tailSurfacesNumber;
}

std::size_t schemi::mesh::innerNumber() const noexcept
{
	return innerSurfacesNumber;
}

std::size_t schemi::mesh::pointNumber() const noexcept
{
	return pointSurfacesNumber;
}

std::size_t schemi::mesh::bottomNumber() const noexcept
{
	return bottomSurfacesNumber;
}

std::size_t schemi::mesh::rightNumber() const noexcept
{
	return rightSurfacesNumber;
}

std::size_t schemi::mesh::leftNumber() const noexcept
{
	return leftSurfacesNumber;
}

std::size_t schemi::mesh::topNumber() const noexcept
{
	return topSurfacesNumber;
}

std::size_t schemi::mesh::nonexistCell() const noexcept
{
	return nonexistentCell;
}

std::size_t schemi::mesh::cellsSize() const noexcept
{
	return cellsNumber;
}

std::size_t schemi::mesh::surfacesSize() const noexcept
{
	return surfacesNumber;
}

schemi::dimensionsEnum schemi::mesh::taskDimension() const noexcept
{
	return taskDim;
}

const std::array<std::size_t, 3>& schemi::mesh::nCells() const noexcept
{
	return n_cells;
}

schemi::scalar schemi::mesh::timestep() const noexcept
{
	return timestep_val;
}

void schemi::mesh::setTimestep(const scalar dt) noexcept
{
	timestep_val = dt;
}

schemi::scalar schemi::mesh::timestepSource() const noexcept
{
	return timestepSource_val;
}

schemi::scalar& schemi::mesh::timestepSourceRef() noexcept
{
	return timestepSource_val;
}

void schemi::mesh::oneDParallelepiped(
		const std::pair<vector, scalar> & vectorOfParallelepiped,
		const std::size_t N_x,
		const std::vector<boundaryConditionType> & commonConditions)
{
	taskDim = dimensionsEnum::task1D;

	n_cells = { N_x, 1, 1 };

	const scalar dx { vectorOfParallelepiped.first.v()[0] / N_x };

	/*Add cells*/
	for (std::size_t i = 0; i < N_x; ++i)
	{
		cellsA.push_back(cubicCell());
		cubicCell & cell = cellsA[i];

		cell.r000_r() = vector(vectorOfParallelepiped.second + dx * i, 0., 0.);
		cell.rX00_r() = vector(vectorOfParallelepiped.second + dx * (i + 1), 0.,
				0.);
		cell.r0Y0_r() = vector(vectorOfParallelepiped.second + dx * i,
				vectorOfParallelepiped.first.v()[1], 0.);
		cell.rXY0_r() = vector(vectorOfParallelepiped.second + dx * (i + 1),
				vectorOfParallelepiped.first.v()[1], 0.);
		cell.r00Z_r() = vector(vectorOfParallelepiped.second + dx * i, 0.,
				vectorOfParallelepiped.first.v()[2]);
		cell.rX0Z_r() = vector(vectorOfParallelepiped.second + dx * (i + 1), 0.,
				vectorOfParallelepiped.first.v()[2]);
		cell.r0YZ_r() = vector(vectorOfParallelepiped.second + dx * i,
				vectorOfParallelepiped.first.v()[1],
				vectorOfParallelepiped.first.v()[2]);
		cell.rXYZ_r() = vector(vectorOfParallelepiped.second + dx * (i + 1),
				vectorOfParallelepiped.first.v()[1],
				vectorOfParallelepiped.first.v()[2]);

		cell.rC_r() = (cell.r000() + cell.rX00() + cell.r0Y0() + cell.rXY0()
				+ cell.r00Z() + cell.rX0Z() + cell.r0YZ() + cell.rXYZ()) / 8.;

		cell.V_r() = (cell.rX00() - cell.r000()).mag()
				* (cell.r0Y0() - cell.r000()).mag()
				* (cell.r00Z() - cell.r000()).mag();
	}

	/*Add tail surface*/
	if (commonConditions[0] == boundaryConditionType::calculated)
	{
		surfacesA.push_back(quadraticSurface());
		surfaceBoundaryCondition.push_back(commonConditions[0]);
		tailSurfacesNumber++;
		quadraticSurface & surface = surfacesA[0];

		surface.r00_r() = vector(vectorOfParallelepiped.second + 0., 0., 0.);
		surface.rX0_r() = vector(vectorOfParallelepiped.second + 0., 0.,
				vectorOfParallelepiped.first.v()[2]);
		surface.r0Y_r() = vector(vectorOfParallelepiped.second + 0.,
				vectorOfParallelepiped.first.v()[1], 0.);
		surface.rXY_r() = vector(vectorOfParallelepiped.second + 0.,
				vectorOfParallelepiped.first.v()[1],
				vectorOfParallelepiped.first.v()[2]);

		surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
				+ surface.rXY()) / 4.;
		surface.S_r() = (surface.rX0() - surface.r00()).mag()
				* (surface.r0Y() - surface.r00()).mag();
	}
	else
		throw exception("Tail surface must be marked <<calculated>>.",
				errorsEnum::meshGenerationError);

	/*Add inner surfaces*/
	for (std::size_t i = 1; i < N_x; ++i)
	{
		surfacesA.push_back(quadraticSurface());
		surfaceBoundaryCondition.push_back(boundaryConditionType::innerSurface);
		innerSurfacesNumber++;
		quadraticSurface & surface = surfacesA[i];

		surface.r00_r() = vector(vectorOfParallelepiped.second + dx * i, 0.,
				0.);
		surface.rX0_r() = vector(vectorOfParallelepiped.second + dx * i, 0.,
				vectorOfParallelepiped.first.v()[2]);
		surface.r0Y_r() = vector(vectorOfParallelepiped.second + dx * i,
				vectorOfParallelepiped.first.v()[1], 0.);
		surface.rXY_r() = vector(vectorOfParallelepiped.second + dx * i,
				vectorOfParallelepiped.first.v()[1],
				vectorOfParallelepiped.first.v()[2]);

		surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
				+ surface.rXY()) / 4.;
		surface.S_r() = (surface.rX0() - surface.r00()).mag()
				* (surface.r0Y() - surface.r00()).mag();
	}

	/*Add point surface*/
	if (commonConditions[1] == boundaryConditionType::calculated)
	{
		surfacesA.push_back(quadraticSurface());
		surfaceBoundaryCondition.push_back(commonConditions[1]);
		pointSurfacesNumber++;
		quadraticSurface & surface = surfacesA[N_x];

		surface.r00_r() = vector(
				vectorOfParallelepiped.second
						+ vectorOfParallelepiped.first.v()[0], 0., 0.);
		surface.rX0_r() = vector(
				vectorOfParallelepiped.second
						+ vectorOfParallelepiped.first.v()[0], 0.,
				vectorOfParallelepiped.first.v()[2]);
		surface.r0Y_r() = vector(
				vectorOfParallelepiped.second
						+ vectorOfParallelepiped.first.v()[0],
				vectorOfParallelepiped.first.v()[1], 0.);
		surface.rXY_r() = vector(
				vectorOfParallelepiped.second
						+ vectorOfParallelepiped.first.v()[0],
				vectorOfParallelepiped.first.v()[1],
				vectorOfParallelepiped.first.v()[2]);

		surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
				+ surface.rXY()) / 4.;
		surface.S_r() = (surface.rX0() - surface.r00()).mag()
				* (surface.r0Y() - surface.r00()).mag();
	}
	else
		throw exception("Point surface must be marked <<calculated>>.",
				errorsEnum::meshGenerationError);

	/*Set surfaces of cell*/
	for (std::size_t i = 0; i < cellsA.size(); ++i)
	{
		surfacesOfCellsA.push_back(std::vector<std::size_t>());
		surfacesOfCellsA[i].resize(2);
		surfacesOfCellsA[i][0] = i;
		surfacesOfCellsA[i][1] = i + 1;
	}

	/*Set environment cells*/
	/*First cell*/
	neighboursOfCellsA.push_back(std::vector<std::size_t>());
	neighboursOfCellsA[0].resize(1);
	neighboursOfCellsA[0][0] = 1;

	/*Other cells*/
	for (std::size_t i = 1; i < (N_x - 1); ++i)
	{
		neighboursOfCellsA.push_back(std::vector<std::size_t>());
		neighboursOfCellsA[i].resize(2);
		neighboursOfCellsA[i][0] = i - 1;
		neighboursOfCellsA[i][1] = i + 1;
	}

	/*Last cell*/
	neighboursOfCellsA.push_back(std::vector<std::size_t>());
	neighboursOfCellsA[N_x - 1].resize(1);
	neighboursOfCellsA[N_x - 1][0] = N_x - 2;

	/*Set surface's owner*/
	{
		/*Tail surface*/
		surfaceOwnerA.push_back(std::size_t());
		surfaceOwnerA[0] = 0;

		/*Inner surfaces*/
		for (std::size_t i = 1; i < N_x; ++i)
		{
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[i] = i - 1;
		}

		/*Point surface*/
		surfaceOwnerA.push_back(std::size_t());
		surfaceOwnerA[N_x] = N_x - 1;
	}

	/*Set surface's neighbour*/
	{
		nonexistentCell = cellsA.size() + 1;

		/*Tail surface*/
		surfaceNeighbourA.push_back(std::size_t());
		surfaceNeighbourA[0] = nonexistentCell;

		/*Inner surfaces*/
		for (std::size_t i = 1; i < N_x; ++i)
		{
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[i] = i;
		}

		/*Point surface*/
		surfaceNeighbourA.push_back(std::size_t());
		surfaceNeighbourA[N_x] = nonexistentCell;
	}

	calculateNormales();

	cellsNumber = cellsA.size();
	surfacesNumber = surfacesA.size();

	initialized = true;
}

void schemi::mesh::twoDParallelepiped(
		const std::pair<vector, scalar> & vectorOfParallelepiped,
		const std::size_t N_x, const std::size_t N_y,
		const std::vector<boundaryConditionType> & commonConditions)
{
	taskDim = dimensionsEnum::task2D;

	n_cells = { N_x, N_y, 1 };

	const scalar dx { vectorOfParallelepiped.first.v()[0] / N_x };
	const scalar dy { vectorOfParallelepiped.first.v()[1] / N_y };

	std::size_t prev;

	/*Add cells*/
	{
		for (std::size_t j = 0; j < N_y; ++j)
			for (std::size_t i = 0; i < N_x; ++i)
			{
				prev = cellsA.size();

				cellsA.push_back(cubicCell());
				cubicCell & cell = cellsA[prev];

				cell.r000_r() = vector(vectorOfParallelepiped.second + dx * i,
						dy * j, 0.);
				cell.rX00_r() = vector(
						vectorOfParallelepiped.second + dx * (i + 1), dy * j,
						0.);
				cell.r0Y0_r() = vector(vectorOfParallelepiped.second + dx * i,
						dy * (j + 1), 0.);
				cell.rXY0_r() = vector(
						vectorOfParallelepiped.second + dx * (i + 1),
						dy * (j + 1), 0.);
				cell.r00Z_r() = vector(vectorOfParallelepiped.second + dx * i,
						dy * j, vectorOfParallelepiped.first.v()[2]);
				cell.rX0Z_r() = vector(
						vectorOfParallelepiped.second + dx * (i + 1), dy * j,
						vectorOfParallelepiped.first.v()[2]);
				cell.r0YZ_r() = vector(vectorOfParallelepiped.second + dx * i,
						dy * (j + 1), vectorOfParallelepiped.first.v()[2]);
				cell.rXYZ_r() = vector(
						vectorOfParallelepiped.second + dx * (i + 1),
						dy * (j + 1), vectorOfParallelepiped.first.v()[2]);

				cell.rC_r() = (cell.r000() + cell.rX00() + cell.r0Y0()
						+ cell.rXY0() + cell.r00Z() + cell.rX0Z() + cell.r0YZ()
						+ cell.rXYZ()) / 8.;

				cell.V_r() = (cell.rX00() - cell.r000()).mag()
						* (cell.r0Y0() - cell.r000()).mag()
						* (cell.r00Z() - cell.r000()).mag();
			}
	}

	/*Add tail surface*/
	{
		if (commonConditions[0] == boundaryConditionType::calculated)
		{
			for (std::size_t j = 0; j < N_y; ++j)
			{
				prev = surfacesA.size();

				surfacesA.push_back(quadraticSurface());
				surfaceBoundaryCondition.push_back(commonConditions[0]);
				tailSurfacesNumber++;
				quadraticSurface & surface = surfacesA[prev];

				surface.r00_r() = vector(vectorOfParallelepiped.second + 0.,
						dy * j, 0.);
				surface.rX0_r() = vector(vectorOfParallelepiped.second + 0.,
						dy * j, vectorOfParallelepiped.first.v()[2]);
				surface.r0Y_r() = vector(vectorOfParallelepiped.second + 0.,
						dy * (j + 1), 0.);
				surface.rXY_r() = vector(vectorOfParallelepiped.second + 0.,
						dy * (j + 1), vectorOfParallelepiped.first.v()[2]);

				surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
						+ surface.rXY()) / 4.;
				surface.S_r() = (surface.rX0() - surface.r00()).mag()
						* (surface.r0Y() - surface.r00()).mag();
			}
		}
		else
			throw exception("Tail surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Add forward inner surfaces*/
	{
		for (std::size_t i = 0; i < N_y * (N_x - 1); ++i)
		{
			const std::size_t layer_n = i / (N_x - 1);
			const std::size_t index_in_layer = (i - layer_n * (N_x - 1));

			if (index_in_layer != (N_x - 1))
			{
				prev = surfacesA.size();

				surfacesA.push_back(quadraticSurface());
				surfaceBoundaryCondition.push_back(
						boundaryConditionType::innerSurface);
				innerSurfacesNumber++;
				quadraticSurface & surface = surfacesA[prev];

				surface.r00_r() = vector(
						vectorOfParallelepiped.second
								+ dx * (index_in_layer + 1), dy * layer_n, 0.);
				surface.rX0_r() = vector(
						vectorOfParallelepiped.second
								+ dx * (index_in_layer + 1), dy * layer_n,
						vectorOfParallelepiped.first.v()[2]);
				surface.r0Y_r() = vector(
						vectorOfParallelepiped.second
								+ dx * (index_in_layer + 1), dy * (layer_n + 1),
						0.);
				surface.rXY_r() = vector(
						vectorOfParallelepiped.second
								+ dx * (index_in_layer + 1), dy * (layer_n + 1),
						vectorOfParallelepiped.first.v()[2]);

				surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
						+ surface.rXY()) / 4.;
				surface.S_r() = (surface.rX0() - surface.r00()).mag()
						* (surface.r0Y() - surface.r00()).mag();
			}
		}
	}

	/*Add left inner surfaces*/
	{
		for (std::size_t i = 0; i < N_x * (N_y - 1); ++i)
		{
			const std::size_t layer_n = i / N_x;
			const std::size_t index_in_layer = (i - layer_n * (N_y - 1));

			prev = surfacesA.size();

			surfacesA.push_back(quadraticSurface());
			surfaceBoundaryCondition.push_back(
					boundaryConditionType::innerSurface);
			innerSurfacesNumber++;
			quadraticSurface & surface = surfacesA[prev];

			surface.r00_r() = vector(
					vectorOfParallelepiped.second + dx * index_in_layer,
					dy * (layer_n + 1), 0.);
			surface.rX0_r() = vector(
					vectorOfParallelepiped.second + dx * (index_in_layer + 1),
					dy * (layer_n + 1), 0.);
			surface.r0Y_r() = vector(
					vectorOfParallelepiped.second + dx * index_in_layer,
					dy * (layer_n + 1), vectorOfParallelepiped.first.v()[2]);
			surface.rXY_r() = vector(
					vectorOfParallelepiped.second + dx * (index_in_layer + 1),
					dy * (layer_n + 1), vectorOfParallelepiped.first.v()[2]);

			surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
					+ surface.rXY()) / 4.;
			surface.S_r() = (surface.rX0() - surface.r00()).mag()
					* (surface.r0Y() - surface.r00()).mag();
		}
	}

	/*Add point surface*/
	{
		if (commonConditions[1] == boundaryConditionType::calculated)
		{
			for (std::size_t j = 0; j < N_y; ++j)
			{
				prev = surfacesA.size();

				surfacesA.push_back(quadraticSurface());
				surfaceBoundaryCondition.push_back(commonConditions[1]);
				pointSurfacesNumber++;
				quadraticSurface & surface = surfacesA[prev];

				surface.r00_r() = vector(
						vectorOfParallelepiped.second
								+ vectorOfParallelepiped.first.v()[0], dy * j,
						0.);
				surface.rX0_r() = vector(
						vectorOfParallelepiped.second
								+ vectorOfParallelepiped.first.v()[0], dy * j,
						vectorOfParallelepiped.first.v()[2]);
				surface.r0Y_r() = vector(
						vectorOfParallelepiped.second
								+ vectorOfParallelepiped.first.v()[0],
						dy * (j + 1), 0.);
				surface.rXY_r() = vector(
						vectorOfParallelepiped.second
								+ vectorOfParallelepiped.first.v()[0],
						dy * (j + 1), vectorOfParallelepiped.first.v()[2]);

				surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
						+ surface.rXY()) / 4.;
				surface.S_r() = (surface.rX0() - surface.r00()).mag()
						* (surface.r0Y() - surface.r00()).mag();
			}
		}
		else
			throw exception("Point surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Add right surface*/
	{
		if (commonConditions[3] == boundaryConditionType::calculated)
		{
			for (std::size_t j = 0; j < N_x; ++j)
			{
				prev = surfacesA.size();

				surfacesA.push_back(quadraticSurface());
				surfaceBoundaryCondition.push_back(commonConditions[3]);
				rightSurfacesNumber++;
				quadraticSurface & surface = surfacesA[prev];

				surface.r00_r() = vector(vectorOfParallelepiped.second + dx * j,
						0., 0.);
				surface.rX0_r() = vector(vectorOfParallelepiped.second + dx * j,
						0., vectorOfParallelepiped.first.v()[2]);
				surface.r0Y_r() = vector(
						vectorOfParallelepiped.second + dx * (j + 1), 0., 0.);
				surface.rXY_r() = vector(
						vectorOfParallelepiped.second + dx * (j + 1), 0.,
						vectorOfParallelepiped.first.v()[2]);

				surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
						+ surface.rXY()) / 4.;
				surface.S_r() = (surface.rX0() - surface.r00()).mag()
						* (surface.r0Y() - surface.r00()).mag();
			}
		}
		else
			throw exception("Right surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Add left surface*/
	{
		if (commonConditions[4] == boundaryConditionType::calculated)
		{
			for (std::size_t j = 0; j < N_x; ++j)
			{
				prev = surfacesA.size();

				surfacesA.push_back(quadraticSurface());
				surfaceBoundaryCondition.push_back(commonConditions[4]);
				leftSurfacesNumber++;
				quadraticSurface & surface = surfacesA[prev];

				surface.r00_r() = vector(vectorOfParallelepiped.second + dx * j,
						vectorOfParallelepiped.first.v()[1], 0.);
				surface.rX0_r() = vector(vectorOfParallelepiped.second + dx * j,
						vectorOfParallelepiped.first.v()[1],
						vectorOfParallelepiped.first.v()[2]);
				surface.r0Y_r() = vector(
						vectorOfParallelepiped.second + dx * (j + 1),
						vectorOfParallelepiped.first.v()[1], 0.);
				surface.rXY_r() = vector(
						vectorOfParallelepiped.second + dx * (j + 1),
						vectorOfParallelepiped.first.v()[1],
						vectorOfParallelepiped.first.v()[2]);

				surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
						+ surface.rXY()) / 4.;
				surface.S_r() = (surface.rX0() - surface.r00()).mag()
						* (surface.r0Y() - surface.r00()).mag();
			}
		}
		else
			throw exception("Left surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Set surfaces of cell*/
	for (std::size_t i = 0; i < cellsA.size(); ++i)
	{
		const std::size_t layer_n = i / N_x;
		const std::size_t index_in_layer = (i - layer_n * N_x);

		surfacesOfCellsA.push_back(std::vector<std::size_t>());
		surfacesOfCellsA[i].resize(4);

		/*First layer*/
		if (layer_n == 0)
		{
			if (index_in_layer == 0)
			{
				surfacesOfCellsA[i][0] = 0; /*tail surface*/
				surfacesOfCellsA[i][1] = tailSurfacesNumber; /*forward inner surface*/
				surfacesOfCellsA[i][2] = tailSurfacesNumber + N_y * (N_x - 1); /*left inner surface*/
				surfacesOfCellsA[i][3] = tailSurfacesNumber
						+ innerSurfacesNumber + pointSurfacesNumber; /*right boundary surface*/
			}
			else if (index_in_layer == (N_x - 1))
			{
				surfacesOfCellsA[i][0] = tailSurfacesNumber + N_x - 2; /*backward (forward) inner surface*/
				surfacesOfCellsA[i][1] = tailSurfacesNumber
						+ innerSurfacesNumber; /*point surface*/
				surfacesOfCellsA[i][2] = tailSurfacesNumber + N_y * (N_x - 1)
						+ N_x - 1; /*left inner surface*/
				surfacesOfCellsA[i][3] = tailSurfacesNumber
						+ innerSurfacesNumber + pointSurfacesNumber + N_x - 1; /*right boundary surface*/
			}
			else
			{
				surfacesOfCellsA[i][0] = tailSurfacesNumber + index_in_layer
						- 1; /*backward (forward) inner surface*/
				surfacesOfCellsA[i][1] = tailSurfacesNumber + index_in_layer; /*forward inner surface*/
				surfacesOfCellsA[i][2] = tailSurfacesNumber + N_y * (N_x - 1)
						+ index_in_layer; /*left inner surface*/
				surfacesOfCellsA[i][3] = tailSurfacesNumber
						+ innerSurfacesNumber + pointSurfacesNumber
						+ index_in_layer; /*right boundary surface*/
			}
		} /*Last layer*/
		else if (layer_n == (N_y - 1))
		{
			if (index_in_layer == 0)
			{
				surfacesOfCellsA[i][0] = tailSurfacesNumber - 1; /*tail surface*/
				surfacesOfCellsA[i][1] = tailSurfacesNumber
						+ (N_y - 1) * (N_x - 1); /*forward inner surface*/
				surfacesOfCellsA[i][2] = tailSurfacesNumber
						+ innerSurfacesNumber + pointSurfacesNumber
						+ rightSurfacesNumber; /*left boundary surface*/
				surfacesOfCellsA[i][3] = tailSurfacesNumber + N_y * (N_x - 1)
						+ (N_y - 2) * N_x; /*right (left) inner surface*/
			}
			else if (index_in_layer == (N_x - 1))
			{
				surfacesOfCellsA[i][0] = tailSurfacesNumber + N_y * (N_x - 1)
						- 1; /*backward (forward) inner surface*/
				surfacesOfCellsA[i][1] = tailSurfacesNumber
						+ innerSurfacesNumber + pointSurfacesNumber - 1; /*point surface*/
				surfacesOfCellsA[i][2] = tailSurfacesNumber
						+ innerSurfacesNumber + pointSurfacesNumber
						+ rightSurfacesNumber + N_x - 1; /*left boundary surface*/
				surfacesOfCellsA[i][3] = tailSurfacesNumber + N_y * (N_x - 1)
						+ (N_y - 1) * N_x - 1; /*right (left) inner surface*/
			}
			else
			{
				surfacesOfCellsA[i][0] = tailSurfacesNumber
						+ (N_y - 1) * (N_x - 1) + index_in_layer - 1; /*backward (forward) inner surface*/
				surfacesOfCellsA[i][1] = tailSurfacesNumber
						+ (N_y - 1) * (N_x - 1) + index_in_layer; /*forward inner surface*/
				surfacesOfCellsA[i][2] = tailSurfacesNumber
						+ innerSurfacesNumber + pointSurfacesNumber
						+ rightSurfacesNumber + index_in_layer; /*left boundary surface*/
				surfacesOfCellsA[i][3] = tailSurfacesNumber + N_y * (N_x - 1)
						+ (N_y - 2) * N_x + index_in_layer; /*right (left) inner surface*/
			}
		} /*Inner layers*/
		else
		{
			if (index_in_layer == 0)
			{
				surfacesOfCellsA[i][0] = layer_n; /*tail surface*/
				surfacesOfCellsA[i][1] = tailSurfacesNumber
						+ layer_n * (N_x - 1); /*forward inner surface*/
				surfacesOfCellsA[i][2] = tailSurfacesNumber + N_y * (N_x - 1)
						+ layer_n * N_x; /*left inner surface*/
				surfacesOfCellsA[i][3] = tailSurfacesNumber + N_y * (N_x - 1)
						+ (layer_n - 1) * N_x; /*right (left) inner surface*/
			}
			else if (index_in_layer == (N_x - 1))
			{
				surfacesOfCellsA[i][0] = tailSurfacesNumber
						+ (layer_n + 1) * (N_x - 1) - 1; /*backward (forward) inner surface*/
				surfacesOfCellsA[i][1] = tailSurfacesNumber
						+ innerSurfacesNumber + layer_n; /*point surface*/
				surfacesOfCellsA[i][2] = tailSurfacesNumber + N_y * (N_x - 1)
						+ (layer_n + 1) * N_x - 1; /*left inner surface*/
				surfacesOfCellsA[i][3] = tailSurfacesNumber + N_y * (N_x - 1)
						+ layer_n * N_x - 1; /*right (left) inner surface*/
			}
			else
			{
				surfacesOfCellsA[i][0] = tailSurfacesNumber
						+ layer_n * (N_x - 1) + index_in_layer - 1; /*backward (forward) inner surface*/
				surfacesOfCellsA[i][1] = tailSurfacesNumber
						+ layer_n * (N_x - 1) + index_in_layer; /*forward inner surface*/
				surfacesOfCellsA[i][2] = tailSurfacesNumber + N_y * (N_x - 1)
						+ layer_n * N_x + index_in_layer; /*left inner surface*/
				surfacesOfCellsA[i][3] = tailSurfacesNumber + N_y * (N_x - 1)
						+ (layer_n - 1) * N_x + index_in_layer; /*right (left) inner surface*/
			}
		}
	}

	/*Set environment cells*/
	for (std::size_t c = 0; c < cellsA.size(); ++c)
	{
		const std::size_t j = c / N_x;
		const std::size_t i = (c - j * N_x);

		prev = neighboursOfCellsA.size();
		neighboursOfCellsA.push_back(std::vector<std::size_t>());

		if (j == 0)
		{
			if (i == 0)
			{
				neighboursOfCellsA[prev].resize(2);
				neighboursOfCellsA[prev][0] = 1;
				neighboursOfCellsA[prev][1] = N_x;
			}
			else if (i == (N_x - 1))
			{
				neighboursOfCellsA[prev].resize(2);
				neighboursOfCellsA[prev][0] = N_x - 2;
				neighboursOfCellsA[prev][1] = 2 * N_x - 1;
			}
			else
			{
				neighboursOfCellsA[prev].resize(3);
				neighboursOfCellsA[prev][0] = i - 1;
				neighboursOfCellsA[prev][1] = i + 1;
				neighboursOfCellsA[prev][2] = i + N_x;
			}
		}
		else if (j == (N_y - 1))
		{
			if (i == 0)
			{
				neighboursOfCellsA[prev].resize(2);
				neighboursOfCellsA[prev][0] = (N_y - 1) * N_x + 1;
				neighboursOfCellsA[prev][1] = (N_y - 2) * N_x;
			}
			else if (i == (N_x - 1))
			{
				neighboursOfCellsA[prev].resize(2);
				neighboursOfCellsA[prev][0] = N_y * N_x - 2;
				neighboursOfCellsA[prev][1] = (N_y - 1) * N_x - 1;
			}
			else
			{
				neighboursOfCellsA[prev].resize(3);
				neighboursOfCellsA[prev][0] = (N_y - 1) * N_x + i - 1;
				neighboursOfCellsA[prev][1] = (N_y - 1) * N_x + i + 1;
				neighboursOfCellsA[prev][2] = (N_y - 2) * N_x + i;
			}
		}
		else
		{
			if (i == 0)
			{
				neighboursOfCellsA[prev].resize(3);
				neighboursOfCellsA[prev][0] = j * N_x + 1;
				neighboursOfCellsA[prev][1] = (j - 1) * N_x;
				neighboursOfCellsA[prev][2] = (j + 1) * N_x;
			}
			else if (i == (N_x - 1))
			{
				neighboursOfCellsA[prev].resize(3);
				neighboursOfCellsA[prev][0] = (j + 1) * N_x - 2;
				neighboursOfCellsA[prev][1] = j * N_x - 1;
				neighboursOfCellsA[prev][2] = (j + 2) * N_x - 1;
			}
			else
			{
				neighboursOfCellsA[prev].resize(4);
				neighboursOfCellsA[prev][0] = j * N_x + i - 1;
				neighboursOfCellsA[prev][1] = j * N_x + i + 1;
				neighboursOfCellsA[prev][2] = (j - 1) * N_x + i;
				neighboursOfCellsA[prev][3] = (j + 1) * N_x + i;
			}
		}
	}

	/*Set surface's owner*/
	{
		/*Tail surfaces*/
		for (std::size_t i = 0; i < N_y; ++i)
		{
			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = i * N_x;
		}

		/*Inner forward surfaces*/
		for (std::size_t c = 0; c < N_y * N_x; ++c)
		{
			const std::size_t j = c / N_x;
			const std::size_t i = c - j * N_x;

			if (i != (N_x - 1))
			{
				prev = surfaceOwnerA.size();
				surfaceOwnerA.push_back(std::size_t());
				surfaceOwnerA[prev] = c;
			}
		}

		/*Inner left surfaces*/
		for (std::size_t c = 0; c < N_y * N_x; ++c)
		{
			const std::size_t j = c / N_x;

			if (j != (N_y - 1))
			{
				prev = surfaceOwnerA.size();
				surfaceOwnerA.push_back(std::size_t());
				surfaceOwnerA[prev] = c;
			}
		}

		/*Point surfaces*/
		for (std::size_t i = 0; i < N_y; ++i)
		{
			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = (i + 1) * N_x - 1;
		}

		/*Right boundary surfaces*/
		for (std::size_t i = 0; i < N_x; ++i)
		{
			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = i;
		}

		/*Left boundary surfaces*/
		for (std::size_t i = 0; i < N_x; ++i)
		{
			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = (N_y - 1) * N_x + i;
		}
	}

	/*Set surface's neighbour*/
	{
		nonexistentCell = cellsA.size() + 1;

		/*Tail surfaces*/
		for (std::size_t i = 0; i < N_y; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}

		/*Inner forward surfaces*/
		for (std::size_t c = 0; c < N_y * N_x; ++c)
		{
			const std::size_t j = c / N_x;
			const std::size_t i = c - j * N_x;

			if (i != 0)
			{
				prev = surfaceNeighbourA.size();
				surfaceNeighbourA.push_back(std::size_t());
				surfaceNeighbourA[prev] = c;
			}
		}

		/*Inner left surfaces*/
		for (std::size_t c = 0; c < N_y * N_x; ++c)
		{
			const std::size_t j = c / N_x;

			if (j != 0)
			{
				prev = surfaceNeighbourA.size();
				surfaceNeighbourA.push_back(std::size_t());
				surfaceNeighbourA[prev] = c;
			}
		}

		/*Point surfaces*/
		for (std::size_t i = 0; i < N_y; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}

		/*Right boundary surfaces*/
		for (std::size_t i = 0; i < N_x; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}

		/*Left boundary surfaces*/
		for (std::size_t i = 0; i < N_x; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}
	}

	calculateNormales();

	cellsNumber = cellsA.size();
	surfacesNumber = surfacesA.size();

	initialized = true;
}

void schemi::mesh::threeDParallelepiped(
		const std::pair<vector, scalar> & vectorOfParallelepiped,
		const std::size_t N_x, const std::size_t N_y, const std::size_t N_z,
		const std::vector<boundaryConditionType> & commonConditions)
{
	taskDim = dimensionsEnum::task3D;

	n_cells = { N_x, N_y, N_z };

	const scalar dx { vectorOfParallelepiped.first.v()[0] / N_x };
	const scalar dy { vectorOfParallelepiped.first.v()[1] / N_y };
	const scalar dz { vectorOfParallelepiped.first.v()[1] / N_z };

	std::size_t prev;

	/*Add cells*/
	{
		for (std::size_t k = 0; k < N_z; ++k)
			for (std::size_t j = 0; j < N_y; ++j)
				for (std::size_t i = 0; i < N_x; ++i)
				{
					prev = cellsA.size();

					cellsA.push_back(cubicCell());
					cubicCell & cell = cellsA[prev];

					cell.r000_r() = vector(
							vectorOfParallelepiped.second + dx * i, dy * j,
							dz * k);
					cell.rX00_r() = vector(
							vectorOfParallelepiped.second + dx * (i + 1),
							dy * j, dz * k);
					cell.r0Y0_r() = vector(
							vectorOfParallelepiped.second + dx * i,
							dy * (j + 1), dz * k);
					cell.rXY0_r() = vector(
							vectorOfParallelepiped.second + dx * (i + 1),
							dy * (j + 1), dz * k);
					cell.r00Z_r() = vector(
							vectorOfParallelepiped.second + dx * i, dy * j,
							dz * (k + 1));
					cell.rX0Z_r() = vector(
							vectorOfParallelepiped.second + dx * (i + 1),
							dy * j, dz * (k + 1));
					cell.r0YZ_r() = vector(
							vectorOfParallelepiped.second + dx * i,
							dy * (j + 1), dz * (k + 1));
					cell.rXYZ_r() = vector(
							vectorOfParallelepiped.second + dx * (i + 1),
							dy * (j + 1), dz * (k + 1));

					cell.rC_r() = (cell.r000() + cell.rX00() + cell.r0Y0()
							+ cell.rXY0() + cell.r00Z() + cell.rX0Z()
							+ cell.r0YZ() + cell.rXYZ()) / 8.;

					cell.V_r() = (cell.rX00() - cell.r000()).mag()
							* (cell.r0Y0() - cell.r000()).mag()
							* (cell.r00Z() - cell.r000()).mag();
				}
	}

	/*Add tail surface*/
	{
		if (commonConditions[0] == boundaryConditionType::calculated)
		{
			for (std::size_t k = 0; k < N_z; ++k)
				for (std::size_t j = 0; j < N_y; ++j)
				{
					prev = surfacesA.size();

					surfacesA.push_back(quadraticSurface());
					surfaceBoundaryCondition.push_back(commonConditions[0]);
					tailSurfacesNumber++;
					quadraticSurface & surface = surfacesA[prev];

					surface.r00_r() = vector(vectorOfParallelepiped.second + 0.,
							dy * j, dz * k);
					surface.rX0_r() = vector(vectorOfParallelepiped.second + 0.,
							dy * j, dz * (k + 1));
					surface.r0Y_r() = vector(vectorOfParallelepiped.second + 0.,
							dy * (j + 1), dz * k);
					surface.rXY_r() = vector(vectorOfParallelepiped.second + 0.,
							dy * (j + 1), dz * (k + 1));

					surface.rC_r() = (surface.r00() + surface.rX0()
							+ surface.r0Y() + surface.rXY()) / 4.;
					surface.S_r() = (surface.rX0() - surface.r00()).mag()
							* (surface.r0Y() - surface.r00()).mag();
				}
		}
		else
			throw exception("Tail surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Add forward inner surfaces*/
	{
		for (std::size_t i = 0; i < N_z * N_y * (N_x - 1); ++i)
		{
			const std::size_t layer_z = i / (N_y * (N_x - 1));
			const std::size_t layer_y = (i - layer_z * N_y * (N_x - 1))
					/ (N_x - 1);
			const std::size_t index_in_layer = (i - layer_z * N_y * (N_x - 1)
					- layer_y * (N_x - 1));

			if (index_in_layer != (N_x - 1))
			{
				prev = surfacesA.size();

				surfacesA.push_back(quadraticSurface());
				surfaceBoundaryCondition.push_back(
						boundaryConditionType::innerSurface);
				innerSurfacesNumber++;
				quadraticSurface & surface = surfacesA[prev];

				surface.r00_r() = vector(
						vectorOfParallelepiped.second
								+ dx * (index_in_layer + 1), dy * layer_y,
						dz * layer_z);
				surface.rX0_r() = vector(
						vectorOfParallelepiped.second
								+ dx * (index_in_layer + 1), dy * layer_y,
						dz * (layer_z + 1));
				surface.r0Y_r() = vector(
						vectorOfParallelepiped.second
								+ dx * (index_in_layer + 1), dy * (layer_y + 1),
						dz * layer_z);
				surface.rXY_r() = vector(
						vectorOfParallelepiped.second
								+ dx * (index_in_layer + 1), dy * (layer_y + 1),
						dz * (layer_z + 1));

				surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
						+ surface.rXY()) / 4.;
				surface.S_r() = (surface.rX0() - surface.r00()).mag()
						* (surface.r0Y() - surface.r00()).mag();
			}
		}
	}

	/*Add left inner surfaces*/
	{
		for (std::size_t i = 0; i < N_z * (N_y - 1) * N_x; ++i)
		{
			const std::size_t layer_z = i / ((N_y - 1) * N_x);
			const std::size_t layer_y = (i - layer_z * (N_y - 1) * N_x) / N_x;
			const std::size_t index_in_layer = (i - layer_z * (N_y - 1) * N_x
					- layer_y * N_x);

			prev = surfacesA.size();

			surfacesA.push_back(quadraticSurface());
			surfaceBoundaryCondition.push_back(
					boundaryConditionType::innerSurface);
			innerSurfacesNumber++;
			quadraticSurface & surface = surfacesA[prev];

			surface.r00_r() = vector(
					vectorOfParallelepiped.second + dx * index_in_layer,
					dy * (layer_y + 1), dz * layer_z);
			surface.rX0_r() = vector(
					vectorOfParallelepiped.second + dx * (index_in_layer + 1),
					dy * (layer_y + 1), dz * layer_z);
			surface.r0Y_r() = vector(
					vectorOfParallelepiped.second + dx * index_in_layer,
					dy * (layer_y + 1), dz * (layer_z + 1));
			surface.rXY_r() = vector(
					vectorOfParallelepiped.second + dx * (index_in_layer + 1),
					dy * (layer_y + 1), dz * (layer_z + 1));

			surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
					+ surface.rXY()) / 4.;
			surface.S_r() = (surface.rX0() - surface.r00()).mag()
					* (surface.r0Y() - surface.r00()).mag();
		}
	}

	/*Add top inner surfaces*/
	{
		for (std::size_t i = 0; i < (N_z - 1) * N_y * N_x; ++i)
		{
			const std::size_t layer_z = i / (N_y * N_x);
			const std::size_t layer_y = (i - layer_z * N_y * N_x) / N_x;
			const std::size_t index_in_layer = (i - layer_z * N_y * N_x
					- layer_y * N_x);

			prev = surfacesA.size();

			surfacesA.push_back(quadraticSurface());
			surfaceBoundaryCondition.push_back(
					boundaryConditionType::innerSurface);
			innerSurfacesNumber++;
			quadraticSurface & surface = surfacesA[prev];

			surface.r00_r() = vector(
					vectorOfParallelepiped.second + dx * index_in_layer,
					dy * layer_y, dz * (layer_z + 1));
			surface.rX0_r() = vector(
					vectorOfParallelepiped.second + dx * (index_in_layer + 1),
					dy * layer_y, dz * (layer_z + 1));
			surface.r0Y_r() = vector(
					vectorOfParallelepiped.second + dx * index_in_layer,
					dy * (layer_y + 1), dz * (layer_z + 1));
			surface.rXY_r() = vector(
					vectorOfParallelepiped.second + dx * (index_in_layer + 1),
					dy * (layer_y + 1), dz * (layer_z + 1));

			surface.rC_r() = (surface.r00() + surface.rX0() + surface.r0Y()
					+ surface.rXY()) / 4.;
			surface.S_r() = (surface.rX0() - surface.r00()).mag()
					* (surface.r0Y() - surface.r00()).mag();
		}
	}

	/*Add point surface*/
	{
		if (commonConditions[1] == boundaryConditionType::calculated)
		{
			for (std::size_t k = 0; k < N_z; ++k)
				for (std::size_t j = 0; j < N_y; ++j)
				{
					prev = surfacesA.size();

					surfacesA.push_back(quadraticSurface());
					surfaceBoundaryCondition.push_back(commonConditions[1]);
					pointSurfacesNumber++;
					quadraticSurface & surface = surfacesA[prev];

					surface.r00_r() = vector(
							vectorOfParallelepiped.second
									+ vectorOfParallelepiped.first.v()[0],
							dy * j, dz * k);
					surface.rX0_r() = vector(
							vectorOfParallelepiped.second
									+ vectorOfParallelepiped.first.v()[0],
							dy * j, dz * (k + 1));
					surface.r0Y_r() = vector(
							vectorOfParallelepiped.second
									+ vectorOfParallelepiped.first.v()[0],
							dy * (j + 1), dz * k);
					surface.rXY_r() = vector(
							vectorOfParallelepiped.second
									+ vectorOfParallelepiped.first.v()[0],
							dy * (j + 1), dz * (k + 1));

					surface.rC_r() = (surface.r00() + surface.rX0()
							+ surface.r0Y() + surface.rXY()) / 4.;
					surface.S_r() = (surface.rX0() - surface.r00()).mag()
							* (surface.r0Y() - surface.r00()).mag();
				}
		}
		else
			throw exception("Point surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Add bottom surface*/
	{
		if (commonConditions[2] == boundaryConditionType::calculated)
		{
			for (std::size_t j = 0; j < N_y; ++j)
				for (std::size_t i = 0; i < N_x; ++i)
				{
					prev = surfacesA.size();

					surfacesA.push_back(quadraticSurface());
					surfaceBoundaryCondition.push_back(commonConditions[2]);
					bottomSurfacesNumber++;
					quadraticSurface & surface = surfacesA[prev];

					surface.r00_r() = vector(
							vectorOfParallelepiped.second + dx * i, dy * j, 0.);
					surface.rX0_r() = vector(
							vectorOfParallelepiped.second + dx * (i + 1),
							dy * j, 0.);
					surface.r0Y_r() = vector(
							vectorOfParallelepiped.second + dx * i,
							dy * (j + 1), 0.);
					surface.rXY_r() = vector(
							vectorOfParallelepiped.second + dx * (i + 1),
							dy * (j + 1), 0.);

					surface.rC_r() = (surface.r00() + surface.rX0()
							+ surface.r0Y() + surface.rXY()) / 4.;
					surface.S_r() = (surface.rX0() - surface.r00()).mag()
							* (surface.r0Y() - surface.r00()).mag();
				}
		}
		else
			throw exception("Bottom surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Add right surface*/
	{
		if (commonConditions[3] == boundaryConditionType::calculated)
		{
			for (std::size_t k = 0; k < N_z; ++k)
				for (std::size_t j = 0; j < N_x; ++j)
				{
					prev = surfacesA.size();

					surfacesA.push_back(quadraticSurface());
					surfaceBoundaryCondition.push_back(commonConditions[3]);
					rightSurfacesNumber++;
					quadraticSurface & surface = surfacesA[prev];

					surface.r00_r() = vector(
							vectorOfParallelepiped.second + dx * j, 0., dz * k);
					surface.rX0_r() = vector(
							vectorOfParallelepiped.second + dx * j, 0.,
							dz * (k + 1));
					surface.r0Y_r() = vector(
							vectorOfParallelepiped.second + dx * (j + 1), 0.,
							dz * k);
					surface.rXY_r() = vector(
							vectorOfParallelepiped.second + dx * (j + 1), 0.,
							dz * (k + 1));

					surface.rC_r() = (surface.r00() + surface.rX0()
							+ surface.r0Y() + surface.rXY()) / 4.;
					surface.S_r() = (surface.rX0() - surface.r00()).mag()
							* (surface.r0Y() - surface.r00()).mag();
				}
		}
		else
			throw exception("Right surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Add left surface*/
	{
		if (commonConditions[4] == boundaryConditionType::calculated)
		{
			for (std::size_t k = 0; k < N_z; ++k)
				for (std::size_t j = 0; j < N_x; ++j)
				{
					prev = surfacesA.size();

					surfacesA.push_back(quadraticSurface());
					surfaceBoundaryCondition.push_back(commonConditions[4]);
					leftSurfacesNumber++;
					quadraticSurface & surface = surfacesA[prev];

					surface.r00_r() = vector(
							vectorOfParallelepiped.second + dx * j,
							vectorOfParallelepiped.first.v()[1], dz * k);
					surface.rX0_r() = vector(
							vectorOfParallelepiped.second + dx * j,
							vectorOfParallelepiped.first.v()[1], dz * (k + 1));
					surface.r0Y_r() = vector(
							vectorOfParallelepiped.second + dx * (j + 1),
							vectorOfParallelepiped.first.v()[1], dz * k);
					surface.rXY_r() = vector(
							vectorOfParallelepiped.second + dx * (j + 1),
							vectorOfParallelepiped.first.v()[1], dz * (k + 1));

					surface.rC_r() = (surface.r00() + surface.rX0()
							+ surface.r0Y() + surface.rXY()) / 4.;
					surface.S_r() = (surface.rX0() - surface.r00()).mag()
							* (surface.r0Y() - surface.r00()).mag();
				}
		}
		else
			throw exception("Left surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Add top surface*/
	{
		if (commonConditions[5] == boundaryConditionType::calculated)
		{
			for (std::size_t j = 0; j < N_y; ++j)
				for (std::size_t i = 0; i < N_x; ++i)
				{
					prev = surfacesA.size();

					surfacesA.push_back(quadraticSurface());
					surfaceBoundaryCondition.push_back(commonConditions[5]);
					topSurfacesNumber++;
					quadraticSurface & surface = surfacesA[prev];

					surface.r00_r() = vector(
							vectorOfParallelepiped.second + dx * i, dy * j,
							vectorOfParallelepiped.first.v()[2]);
					surface.rX0_r() = vector(
							vectorOfParallelepiped.second + dx * (i + 1),
							dy * j, vectorOfParallelepiped.first.v()[2]);
					surface.r0Y_r() = vector(
							vectorOfParallelepiped.second + dx * i,
							dy * (j + 1), vectorOfParallelepiped.first.v()[2]);
					surface.rXY_r() = vector(
							vectorOfParallelepiped.second + dx * (i + 1),
							dy * (j + 1), vectorOfParallelepiped.first.v()[2]);

					surface.rC_r() = (surface.r00() + surface.rX0()
							+ surface.r0Y() + surface.rXY()) / 4.;
					surface.S_r() = (surface.rX0() - surface.r00()).mag()
							* (surface.r0Y() - surface.r00()).mag();
				}
		}
		else
			throw exception("Top surface must be marked <<calculated>>.",
					errorsEnum::meshGenerationError);
	}

	/*Set surfaces of cell*/
	for (std::size_t i = 0; i < cellsA.size(); ++i)
	{
		const std::size_t layer_z = i / (N_y * N_x);
		const std::size_t layer_y = (i - layer_z * N_y * N_x) / N_x;
		const std::size_t index_in_layer = (i - layer_z * N_y * N_x
				- layer_y * N_x);

		surfacesOfCellsA.push_back(std::vector<std::size_t>());
		surfacesOfCellsA[i].resize(6);

		if (layer_z == 0)
		{
			if (layer_y == 0)
			{
				if (index_in_layer == 0)
				{
					surfacesOfCellsA[i][0] = 0; /*tail surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1); /*forward inner surfaces*/
					/*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber; /*right boundary surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber; /*bottom boundary surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x; /*left inner surfaces*/
					/*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber + N_x - 2; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ innerSurfacesNumber; /*point surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_x - 1; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + N_x - 1; /*right boundary surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber + N_x
							- 1; /*bottom boundary surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ N_x - 1; /*top inner surface*/
				}
				else
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber + index_in_layer
							- 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ index_in_layer; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ index_in_layer; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + index_in_layer; /*right boundary surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ index_in_layer; /*bottom boundary surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ index_in_layer; /*top inner surface*/
				}
			}
			else if (layer_y == (N_y - 1))
			{
				if (index_in_layer == 0)
				{
					surfacesOfCellsA[i][0] = N_y - 1; /*tails surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ (N_y - 1) * (N_x - 1); /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber; /*left boundary surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_y - 2) * N_x; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ (N_y - 1) * N_x; /*bottom boundary surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_y - 1) * N_x; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ (N_y - 1) * (N_x - 1) + N_x - 2; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ innerSurfacesNumber + N_y - 1; /*point surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber + N_x
							- 1; /*left boundary surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_y - 2) * N_x + N_x - 1; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ (N_y - 1) * N_x + N_x - 1; /*bottom boundary surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_y - 1) * N_x + N_x - 1; /*top inner surface*/
				}
				else
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ (N_y - 1) * (N_x - 1) + index_in_layer - 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ (N_y - 1) * (N_x - 1) + index_in_layer; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ index_in_layer; /*left boundary surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_y - 2) * N_x + index_in_layer; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ (N_y - 1) * N_x + index_in_layer; /*bottom boundary surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_y - 1) * N_x + index_in_layer; /*top inner surface*/
				}
			}
			else
			{
				if (index_in_layer == 0)
				{
					surfacesOfCellsA[i][0] = layer_y; /*tail surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ layer_y * (N_x - 1); /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_y * N_x; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (layer_y - 1) * N_x; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ layer_y * N_x; /*bottom boundary surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_y * N_x; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ layer_y * (N_x - 1) + N_x - 2; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ innerSurfacesNumber + layer_y; /*point surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_y * N_x + N_x - 1; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (layer_y - 1) * N_x + N_x - 1; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ layer_y * N_x + N_x - 1; /*bottom boundary surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_y * N_x + N_x - 1; /*top inner surface*/
				}
				else
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ layer_y * (N_x - 1) + index_in_layer - 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ layer_y * (N_x - 1) + index_in_layer; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_y * N_x + index_in_layer; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ layer_y * N_x + index_in_layer; /*bottom boundary surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
			}
		}
		else if (layer_z == (N_z - 1))
		{
			if (layer_y == 0)
			{
				if (index_in_layer == 0)
				{
					surfacesOfCellsA[i][0] = (N_z - 1) * N_y; /*tail surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1); /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + (N_z - 1) * N_x; /*right boundary surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_z - 2) * N_y * N_x; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ leftSurfacesNumber; /*top boundary surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1) + N_x - 1 - 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ innerSurfacesNumber + (N_z - 1) * N_y; /*point surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + N_x - 1; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + (N_z - 1) * N_x + N_x - 1; /*right boundary surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_z - 2) * N_y * N_x + N_x - 1; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ leftSurfacesNumber + N_x - 1; /*top boundary surface*/
				}
				else
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1) + index_in_layer - 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1) + index_in_layer; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + index_in_layer; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + (N_z - 1) * N_x
							+ index_in_layer; /*right boundary surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_z - 2) * N_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ leftSurfacesNumber + index_in_layer; /*top boundary surface*/
				}
			}
			else if (layer_y == (N_y - 1))
			{
				if (index_in_layer == 0)
				{
					surfacesOfCellsA[i][0] = (N_z - 1) * N_y + N_y - 1; /*tail surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1)
							+ (N_y - 1) * (N_x - 1); /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ (N_z - 1) * N_x; /*left boundary surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + (N_y - 2) * N_x; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_z - 2) * N_y * N_x + (N_y - 1) * N_x; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ leftSurfacesNumber + (N_y - 1) * N_x; /*top boundary surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1)
							+ (N_y - 1) * (N_x - 1) + N_x - 2; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ innerSurfacesNumber + (N_z - 1) * N_y + N_y - 1; /*point surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ (N_z - 1) * N_x + N_x - 1; /*left boundary surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + (N_y - 2) * N_x
							+ N_x - 1; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_z - 2) * N_y * N_x + (N_y - 1) * N_x + N_x - 1; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ leftSurfacesNumber + (N_y - 1) * N_x + N_x - 1; /*top boundary surface*/
				}
				else
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1)
							+ (N_y - 1) * (N_x - 1) + index_in_layer - 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1)
							+ (N_y - 1) * (N_x - 1) + index_in_layer; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ (N_z - 1) * N_x + index_in_layer; /*left boundary surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + (N_y - 2) * N_x
							+ index_in_layer; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_z - 2) * N_y * N_x + (N_y - 1) * N_x
							+ index_in_layer; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ leftSurfacesNumber + (N_y - 1) * N_x
							+ index_in_layer; /*top boundary surface*/
				}
			}
			else
			{
				if (index_in_layer == 0)
				{
					surfacesOfCellsA[i][0] = (N_z - 1) * N_y + layer_y; /*tail surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1) + layer_y * (N_x - 1); /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + layer_y * N_x; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + (layer_y - 1) * N_x; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_z - 2) * N_y * N_x + layer_y * N_x; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ leftSurfacesNumber + layer_y * N_x; /*top boundary surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1) + layer_y * (N_x - 1)
							+ N_x - 2; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ innerSurfacesNumber + (N_z - 1) * N_y + layer_y; /*point surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + layer_y * N_x + N_x
							- 1; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + (layer_y - 1) * N_x
							+ N_x - 1; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_z - 2) * N_y * N_x + layer_y * N_x + N_x - 1; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ leftSurfacesNumber + layer_y * N_x + N_x - 1; /*top boundary surface*/
				}
				else
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1) + layer_y * (N_x - 1)
							+ index_in_layer - 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ (N_z - 1) * N_y * (N_x - 1) + layer_y * (N_x - 1)
							+ index_in_layer; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + layer_y * N_x
							+ index_in_layer; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ (N_z - 1) * (N_y - 1) * N_x + (layer_y - 1) * N_x
							+ index_in_layer; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (N_z - 2) * N_y * N_x + layer_y * N_x
							+ index_in_layer; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ leftSurfacesNumber + layer_y * N_x
							+ index_in_layer; /*top boundary surface*/
				}
			}
		}
		else
		{
			if (layer_y == 0)
			{
				if (index_in_layer == 0)
				{
					surfacesOfCellsA[i][0] = layer_z * N_y; /*tail surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1); /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + layer_z * N_x; /*right boundary surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (layer_z - 1) * N_y * N_x; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_z * N_y * N_x; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + N_x - 2; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ innerSurfacesNumber + layer_z * N_y; /*point surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + N_x - 1; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + layer_z * N_x + N_x - 1; /*right boundary surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (layer_z - 1) * N_y * N_x + N_x - 1; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_z * N_y * N_x + N_x - 1; /*top inner surface*/
				}
				else
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + index_in_layer - 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + index_in_layer; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + index_in_layer; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + layer_z * N_x
							+ index_in_layer; /*right boundary surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (layer_z - 1) * N_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_z * N_y * N_x + index_in_layer; /*top inner surface*/
				}
			}
			else if (layer_y == (N_y - 1))
			{
				if (index_in_layer == 0)
				{
					surfacesOfCellsA[i][0] = layer_z * N_y + N_y - 1; /*tail surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + (N_y - 1) * (N_x - 1); /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ layer_z * N_x; /*left boundary surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + (N_y - 2) * N_x; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (layer_z - 1) * N_y * N_x + (N_y - 1) * N_x; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_z * N_y * N_x + (N_y - 1) * N_x; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + (N_y - 1) * (N_x - 1)
							+ N_x - 2; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ innerSurfacesNumber + layer_z * N_y + N_y - 1; /*point surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ layer_z * N_x + N_x - 1; /*left boundary surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + (N_y - 2) * N_x + N_x
							- 1; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (layer_z - 1) * N_y * N_x + (N_y - 1) * N_x + N_x
							- 1; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_z * N_y * N_x + (N_y - 1) * N_x + N_x - 1; /*top inner surface*/
				}
				else
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + (N_y - 1) * (N_x - 1)
							+ index_in_layer - 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + (N_y - 1) * (N_x - 1)
							+ index_in_layer; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ innerSurfacesNumber + pointSurfacesNumber
							+ bottomSurfacesNumber + rightSurfacesNumber
							+ layer_z * N_x + index_in_layer; /*left boundary surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + (N_y - 2) * N_x
							+ index_in_layer; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (layer_z - 1) * N_y * N_x + (N_y - 1) * N_x
							+ index_in_layer; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_z * N_y * N_x + (N_y - 1) * N_x
							+ index_in_layer; /*top inner surface*/
				}
			}
			else
			{
				if (index_in_layer == 0)
				{
					surfacesOfCellsA[i][0] = layer_z * N_y + layer_y; /*tail surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + layer_y * (N_x - 1); /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + layer_y * N_x; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + (layer_y - 1) * N_x; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (layer_z - 1) * N_y * N_x + layer_y * N_x; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_z * N_y * N_x + layer_y * N_x; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + layer_y * (N_x - 1)
							+ N_x - 2; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ innerSurfacesNumber + layer_z * N_y + layer_y; /*point surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + layer_y * N_x + N_x
							- 1; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + (layer_y - 1) * N_x
							+ N_x - 1; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (layer_z - 1) * N_y * N_x + layer_y * N_x + N_x
							- 1; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_z * N_y * N_x + layer_y * N_x + N_x - 1; /*top inner surface*/
				}
				else
				{
					surfacesOfCellsA[i][0] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + layer_y * (N_x - 1)
							+ index_in_layer - 1; /*backward (forward) inner surface*/
					surfacesOfCellsA[i][1] = tailSurfacesNumber
							+ layer_z * N_y * (N_x - 1) + layer_y * (N_x - 1)
							+ index_in_layer; /*forward inner surface*/
					surfacesOfCellsA[i][2] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + layer_y * N_x
							+ index_in_layer; /*left inner surface*/
					surfacesOfCellsA[i][3] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ layer_z * (N_y - 1) * N_x + (layer_y - 1) * N_x
							+ index_in_layer; /*right (left) inner surface*/
					surfacesOfCellsA[i][4] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ (layer_z - 1) * N_y * N_x + layer_y * N_x
							+ index_in_layer; /*bottom (top) inner surface*/
					surfacesOfCellsA[i][5] = tailSurfacesNumber
							+ N_z * N_y * (N_x - 1) /*forward inner surfaces*/
							+ N_z * (N_y - 1) * N_x /*left inner surfaces*/
							+ layer_z * N_y * N_x + layer_y * N_x
							+ index_in_layer; /*top inner surface*/
				}
			}
		}
	}

	/*Set environment cells*/
	for (std::size_t c = 0; c < cellsA.size(); ++c)
	{
		const std::size_t layer_z = c / (N_y * N_x);
		const std::size_t layer_y = (c - layer_z * N_y * N_x) / N_x;
		const std::size_t index_in_layer = (c - layer_z * N_y * N_x
				- layer_y * N_x);

		prev = neighboursOfCellsA.size();
		neighboursOfCellsA.push_back(std::vector<std::size_t>());

		if (layer_z == 0)
		{
			if (layer_y == 0)
			{
				if (index_in_layer == 0)
				{
					neighboursOfCellsA[prev].resize(3);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					neighboursOfCellsA[prev].resize(3);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
			}
			else if (layer_y == (N_y - 1))
			{
				if (index_in_layer == 0)
				{
					neighboursOfCellsA[prev].resize(3);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					neighboursOfCellsA[prev].resize(3);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
			}
			else
			{
				if (index_in_layer == 0)
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else
				{
					neighboursOfCellsA[prev].resize(5);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][3] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][4] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
			}
		}
		else if (layer_z == (N_z - 1))
		{
			if (layer_y == 0)
			{
				if (index_in_layer == 0)
				{
					neighboursOfCellsA[prev].resize(3);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					neighboursOfCellsA[prev].resize(3);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
				}
				else
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
				}
			}
			else if (layer_y == (N_y - 1))
			{
				if (index_in_layer == 0)
				{
					neighboursOfCellsA[prev].resize(3);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					neighboursOfCellsA[prev].resize(3);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
				}
				else
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
				}

			}
			else
			{
				if (index_in_layer == 0)
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
				}
				else
				{
					neighboursOfCellsA[prev].resize(5);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][3] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][4] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
				}
			}
		}
		else
		{
			if (layer_y == 0)
			{
				if (index_in_layer == 0)
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else
				{
					neighboursOfCellsA[prev].resize(5);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					neighboursOfCellsA[prev][4] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
			}
			else if (layer_y == (N_y - 1))
			{
				if (index_in_layer == 0)
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					neighboursOfCellsA[prev].resize(4);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][2] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else
				{
					neighboursOfCellsA[prev].resize(5);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					neighboursOfCellsA[prev][4] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
			}
			else
			{
				if (index_in_layer == 0)
				{
					neighboursOfCellsA[prev].resize(5);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					neighboursOfCellsA[prev][4] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else if (index_in_layer == (N_x - 1))
				{
					neighboursOfCellsA[prev].resize(5);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][3] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					neighboursOfCellsA[prev][4] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
				else
				{
					neighboursOfCellsA[prev].resize(6);
					neighboursOfCellsA[prev][0] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer - 1; /*backward (forward) inner surface*/
					neighboursOfCellsA[prev][1] = layer_z * N_y * N_x
							+ layer_y * N_x + index_in_layer + 1; /*forward inner surface*/
					neighboursOfCellsA[prev][2] = layer_z * N_y * N_x
							+ (layer_y + 1) * N_x + index_in_layer; /*left inner surface*/
					neighboursOfCellsA[prev][3] = layer_z * N_y * N_x
							+ (layer_y - 1) * N_x + index_in_layer; /*right (left) inner surface*/
					neighboursOfCellsA[prev][4] = (layer_z - 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*bottom (top) inner surface*/
					neighboursOfCellsA[prev][5] = (layer_z + 1) * N_y * N_x
							+ layer_y * N_x + index_in_layer; /*top inner surface*/
				}
			}
		}
	}

	/*Set surface's owner*/
	{
		/*Tail surfaces*/
		for (std::size_t i = 0; i < N_z * N_y; ++i)
		{
			const std::size_t layer_z = i / N_y;
			const std::size_t layer_y = i - layer_z * N_y;

			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = layer_z * N_y * N_x + layer_y * N_x;
		}

		/*Inner forward surfaces*/
		for (std::size_t c = 0; c < N_z * N_y * N_x; ++c)
		{
			const std::size_t layer_z = c / (N_y * N_x);
			const std::size_t layer_y = (c - layer_z * N_y * N_x) / N_x;
			const std::size_t index_in_layer = (c - layer_z * N_y * N_x
					- layer_y * N_x);

			if (index_in_layer != (N_x - 1))
			{
				prev = surfaceOwnerA.size();
				surfaceOwnerA.push_back(std::size_t());
				surfaceOwnerA[prev] = c;
			}
		}

		/*Inner left surfaces*/
		for (std::size_t c = 0; c < N_z * N_y * N_x; ++c)
		{
			const std::size_t layer_z = c / (N_y * N_x);
			const std::size_t layer_y = (c - layer_z * N_y * N_x) / N_x;

			if (layer_y != (N_y - 1))
			{
				prev = surfaceOwnerA.size();
				surfaceOwnerA.push_back(std::size_t());
				surfaceOwnerA[prev] = c;
			}
		}

		/*Inner top surfaces*/
		for (std::size_t c = 0; c < N_z * N_y * N_x; ++c)
		{
			const std::size_t layer_z = c / (N_y * N_x);

			if (layer_z != (N_z - 1))
			{
				prev = surfaceOwnerA.size();
				surfaceOwnerA.push_back(std::size_t());
				surfaceOwnerA[prev] = c;
			}
		}

		/*Point surfaces*/
		for (std::size_t i = 0; i < N_z * N_y; ++i)
		{
			const std::size_t layer_z = i / N_y;
			const std::size_t layer_y = i - layer_z * N_y;

			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = layer_z * N_y * N_x + (layer_y + 1) * N_x - 1;
		}

		/*Bottom boundary surfaces*/
		for (std::size_t i = 0; i < N_y * N_x; ++i)
		{
			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = i;
		}

		/*Right boundary surfaces*/
		for (std::size_t i = 0; i < N_z * N_x; ++i)
		{
			const std::size_t layer_z = i / N_x;
			const std::size_t index_in_layer = i - layer_z * N_x;

			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = layer_z * N_y * N_x + index_in_layer;
		}

		/*Left boundary surfaces*/
		for (std::size_t i = 0; i < N_z * N_x; ++i)
		{
			const std::size_t layer_z = i / N_x;
			const std::size_t index_in_layer = i - layer_z * N_x;

			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = layer_z * N_y * N_x + (N_y - 1) * N_x
					+ index_in_layer;
		}

		/*Top boundary surfaces*/
		for (std::size_t i = 0; i < N_y * N_x; ++i)
		{
			prev = surfaceOwnerA.size();
			surfaceOwnerA.push_back(std::size_t());
			surfaceOwnerA[prev] = i + (N_z - 1) * N_y * N_x;
		}
	}

	/*Set surface's neighbour*/
	{
		nonexistentCell = cellsA.size() + 1;

		/*Tail surfaces*/
		for (std::size_t i = 0; i < N_z * N_y; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}

		/*Inner forward surfaces*/
		for (std::size_t c = 0; c < N_z * N_y * N_x; ++c)
		{
			const std::size_t layer_z = c / (N_y * N_x);
			const std::size_t layer_y = (c - layer_z * N_y * N_x) / N_x;
			const std::size_t index_in_layer = (c - layer_z * N_y * N_x
					- layer_y * N_x);

			if (index_in_layer != 0)
			{
				prev = surfaceNeighbourA.size();
				surfaceNeighbourA.push_back(std::size_t());
				surfaceNeighbourA[prev] = c;
			}
		}

		/*Inner left surfaces*/
		for (std::size_t c = 0; c < N_z * N_y * N_x; ++c)
		{
			const std::size_t layer_z = c / (N_y * N_x);
			const std::size_t layer_y = (c - layer_z * N_y * N_x) / N_x;

			if (layer_y != 0)
			{
				prev = surfaceNeighbourA.size();
				surfaceNeighbourA.push_back(std::size_t());
				surfaceNeighbourA[prev] = c;
			}
		}

		/*Inner top surfaces*/
		for (std::size_t c = 0; c < N_z * N_y * N_x; ++c)
		{
			const std::size_t layer_z = c / (N_y * N_x);

			if (layer_z != 0)
			{
				prev = surfaceNeighbourA.size();
				surfaceNeighbourA.push_back(std::size_t());
				surfaceNeighbourA[prev] = c;
			}
		}

		/*Point surfaces*/
		for (std::size_t i = 0; i < N_z * N_y; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}

		/*Bottom boundary surfaces*/
		for (std::size_t i = 0; i < N_y * N_x; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}

		/*Right boundary surfaces*/
		for (std::size_t i = 0; i < N_z * N_x; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}

		/*Left boundary surfaces*/
		for (std::size_t i = 0; i < N_z * N_x; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}

		/*Top boundary surfaces*/
		for (std::size_t i = 0; i < N_y * N_x; ++i)
		{
			prev = surfaceNeighbourA.size();
			surfaceNeighbourA.push_back(std::size_t());
			surfaceNeighbourA[prev] = nonexistentCell;
		}
	}

	calculateNormales();

	cellsNumber = cellsA.size();
	surfacesNumber = surfacesA.size();

	initialized = true;
}

schemi::mesh * schemi::mesh::pInstance = nullptr;
