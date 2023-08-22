/*
 * TVDLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "TVDLimiter.hpp"

schemi::volumeField<schemi::vector> schemi::TVDLimiter(
		const volumeField<vector> & gradient, const volumeField<scalar> & value,
		const abstractLimiter & limiterObjectP,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	auto & mesh { gradient.meshRef() };

	volumeField<vector> retField { mesh, vector(0) };

	for (std::size_t i = 0; i < gradient.size(); ++i)
	{
		vector r;

		const vector & cellVector(mesh.cells()[i].rC());

		std::vector<vector> neighbourGradientsOfSurfaces(
				mesh.neighboursOfCells()[i].size());

		for (std::size_t cx = 0; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			const std::size_t cellIndex { mesh.neighboursOfCells()[i][cx] };

			const vector deltaVector(mesh.cells()[cellIndex].rC() - cellVector);

			const scalar deltaVectorMag { deltaVector.mag() };

			const vector deltaVectorNorm(deltaVector / deltaVectorMag);

			for (std::size_t j = 0; j < vector::vsize; ++j)
			{
				neighbourGradientsOfSurfaces[cx].v_r()[j] =
						(value.ref()[cellIndex] - value.ref()[i])
								/ deltaVectorMag * deltaVectorNorm.v()[j];
			}
		}

		if (mesh.neighboursOfCells()[i].size()
				< mesh.surfacesOfCells()[i].size())
		{
			for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
			{
				const std::size_t surfaceIndex { mesh.surfacesOfCells()[i][j] };

				if (value.boundCond()[surfaceIndex].first
						!= boundaryConditionType::innerSurface)
				{
					const vector deltaVector(
							(mesh.surfaces()[surfaceIndex].rC() - cellVector)
									* 2);

					const scalar deltaVectorMag { deltaVector.mag() };

					const vector deltaVectorNorm(deltaVector / deltaVectorMag);

					const scalar boundaryCellValue =
							bncCalc.boundaryConditionValueCell(value.ref()[i],
									value.boundCond()[surfaceIndex],

									i, surfaceIndex, compt);

					neighbourGradientsOfSurfaces.push_back(
							vector(
									(boundaryCellValue - value.ref()[i])
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i])
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i])
											/ deltaVectorMag
											* deltaVectorNorm.v()[2]));
				}
			}
		}
		else if (mesh.neighboursOfCells()[i].size()
				> mesh.surfacesOfCells()[i].size())
			throw exception(
					"Number of surfaces somehow lesser than number of neighbour cells.",
					errors::systemError);

		vector min(neighbourGradientsOfSurfaces[0]), max(
				neighbourGradientsOfSurfaces[0]);

		for (std::size_t cx = 1; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			for (std::size_t j = 0; j < vector::vsize; ++j)
				min.v_r()[j] = std::min(min.v()[j],
						neighbourGradientsOfSurfaces[cx].v()[j]);

			for (std::size_t j = 0; j < vector::vsize; ++j)
				max.v_r()[j] = std::max(max.v()[j],
						neighbourGradientsOfSurfaces[cx].v()[j]);
		}

		for (std::size_t j = 0; j < vector::vsize; ++j)
			r.v_r()[j] = min.v()[j] / (max.v()[j] + stabilizator);

		retField.ref_r()[i] = limiterObjectP.calculate(r, gradient.ref()[i]);
	}

	return retField;
}

schemi::volumeField<schemi::tensor> schemi::TVDLimiter(
		const volumeField<tensor> & gradient, const volumeField<vector> & value,
		const abstractLimiter & limiterObjectP,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	auto & mesh { gradient.meshRef() };

	volumeField<tensor> retField(mesh, tensor(0));

	for (std::size_t i = 0; i < gradient.size(); ++i)
	{
		tensor r;

		const vector & cellVector(mesh.cells()[i].rC());

		std::vector<tensor> neighbourGradientsOfSurfaces(
				mesh.neighboursOfCells()[i].size());

		for (std::size_t cx = 0; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			const std::size_t cellIndex { mesh.neighboursOfCells()[i][cx] };

			const vector deltaVector(mesh.cells()[cellIndex].rC() - cellVector);

			const scalar deltaVectorMag { deltaVector.mag() };

			const vector deltaVectorNorm(deltaVector / deltaVectorMag);

			for (std::size_t j = 0; j < vector::vsize; ++j)
			{
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[0]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 1] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[1]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 2] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[2]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
			}
		}

		if (mesh.neighboursOfCells()[i].size()
				< mesh.surfacesOfCells()[i].size())
		{
			for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
			{
				const std::size_t surfaceIndex = mesh.surfacesOfCells()[i][j];

				if (value.boundCond()[surfaceIndex].first
						!= boundaryConditionType::innerSurface)
				{
					const vector deltaVector(
							(mesh.surfaces()[surfaceIndex].rC() - cellVector)
									* 2);

					const scalar deltaVectorMag { deltaVector.mag() };

					const vector deltaVectorNorm(deltaVector / deltaVectorMag);

					const vector boundaryCellValue(
							bncCalc.boundaryConditionValueCell(value.ref()[i],
									value.boundCond()[surfaceIndex], i,
									surfaceIndex, compt));

					neighbourGradientsOfSurfaces.push_back(
							tensor(
									(boundaryCellValue - value.ref()[i]).v()[0]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[1]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[2]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[0]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[1]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[2]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[0]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[1]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[2]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2]));
				}
			}
		}
		else if (mesh.neighboursOfCells()[i].size()
				> mesh.surfacesOfCells()[i].size())
			throw exception(
					"Number of surfaces somehow lesser than number of neighbour cells.",
					errors::systemError);

		tensor min(neighbourGradientsOfSurfaces[0]), max(
				neighbourGradientsOfSurfaces[0]);

		for (std::size_t cx = 1; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			for (std::size_t j = 0; j < tensor::vsize; ++j)
				min.v_r()[j] = std::min(min.v()[j],
						neighbourGradientsOfSurfaces[cx].v()[j]);

			for (std::size_t j = 0; j < tensor::vsize; ++j)
				max.v_r()[j] = std::max(max.v()[j],
						neighbourGradientsOfSurfaces[cx].v()[j]);
		}

		for (std::size_t j = 0; j < tensor::vsize; ++j)
			r.v_r()[j] = min.v()[j] / (max.v()[j] + stabilizator);

		retField.ref_r()[i] = limiterObjectP.calculate(r, gradient.ref()[i]);
	}

	return retField;
}

schemi::volumeField<schemi::tensor3> schemi::TVDLimiter(
		const volumeField<tensor3> & gradient,
		const volumeField<tensor> & value,
		const abstractLimiter & limiterObjectP,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	auto & mesh { gradient.meshRef() };

	volumeField<tensor3> retField(mesh, tensor3(0));

	for (std::size_t i = 0; i < gradient.size(); ++i)
	{
		tensor3 r;

		const vector & cellVector(mesh.cells()[i].rC());

		std::vector<tensor3> neighbourGradientsOfSurfaces(
				mesh.neighboursOfCells()[i].size());

		for (std::size_t cx = 0; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			const std::size_t cellIndex { mesh.neighboursOfCells()[i][cx] };

			const vector deltaVector(mesh.cells()[cellIndex].rC() - cellVector);

			const scalar deltaVectorMag { deltaVector.mag() };

			const vector deltaVectorNorm(deltaVector / deltaVectorMag);

			for (std::size_t j = 0; j < vector::vsize; ++j)
			{
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[0]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 1] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[1]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 2] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[2]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 3] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[3]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 4] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[4]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 5] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[5]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 6] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[6]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 7] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[7]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
				neighbourGradientsOfSurfaces[cx].v_r()[3 * j + 8] =
						(value.ref()[cellIndex] - value.ref()[i]).v()[8]
								/ deltaVectorMag * deltaVectorNorm.v()[j];
			}
		}

		if (mesh.neighboursOfCells()[i].size()
				< mesh.surfacesOfCells()[i].size())
		{
			for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
			{
				const std::size_t surfaceIndex = mesh.surfacesOfCells()[i][j];

				if (value.boundCond()[surfaceIndex].first
						!= boundaryConditionType::innerSurface)
				{
					const vector deltaVector(
							(mesh.surfaces()[surfaceIndex].rC() - cellVector)
									* 2);

					const scalar deltaVectorMag { deltaVector.mag() };

					const vector deltaVectorNorm(deltaVector / deltaVectorMag);

					const tensor boundaryCellValue(
							bncCalc.boundaryConditionValueCell(value.ref()[i],
									value.boundCond()[surfaceIndex], i,
									surfaceIndex, compt));

					neighbourGradientsOfSurfaces.push_back(
							tensor3(
									(boundaryCellValue - value.ref()[i]).v()[0]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[1]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[2]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[3]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[4]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[5]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[6]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[7]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],
									(boundaryCellValue - value.ref()[i]).v()[8]
											/ deltaVectorMag
											* deltaVectorNorm.v()[0],

									(boundaryCellValue - value.ref()[i]).v()[0]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[1]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[2]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[3]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[4]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[5]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[6]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[7]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],
									(boundaryCellValue - value.ref()[i]).v()[8]
											/ deltaVectorMag
											* deltaVectorNorm.v()[1],

									(boundaryCellValue - value.ref()[i]).v()[0]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[1]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[2]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[3]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[4]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[5]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[6]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[7]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2],
									(boundaryCellValue - value.ref()[i]).v()[8]
											/ deltaVectorMag
											* deltaVectorNorm.v()[2]));
				}
			}
		}
		else if (mesh.neighboursOfCells()[i].size()
				> mesh.surfacesOfCells()[i].size())
			throw exception(
					"Number of surfaces somehow lesser than number of neighbour cells.",
					errors::systemError);

		tensor3 min(neighbourGradientsOfSurfaces[0]), max(
				neighbourGradientsOfSurfaces[0]);

		for (std::size_t cx = 1; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			for (std::size_t j = 0; j < tensor3::vsize; ++j)
				min.v_r()[j] = std::min(min.v()[j],
						neighbourGradientsOfSurfaces[cx].v()[j]);

			for (std::size_t j = 0; j < tensor3::vsize; ++j)
				max.v_r()[j] = std::max(max.v()[j],
						neighbourGradientsOfSurfaces[cx].v()[j]);
		}

		for (std::size_t j = 0; j < tensor3::vsize; ++j)
			r.v_r()[j] = min.v()[j] / (max.v()[j] + stabilizator);

		retField.ref_r()[i] = limiterObjectP.calculate(r, gradient.ref()[i]);
	}

	return retField;
}
