/*
 * TVDLimiter.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "TVDLimiter.hpp"

#include <algorithm>

schemi::volumeField<schemi::vector> schemi::TVDLimiter(
		const volumeField<vector> & gradient, const volumeField<scalar> & value,
		const abstractLimiter & limiterObjectP,
		const boundaryConditionValue & bncCalc, const MPIHandler & parallel,
		const std::size_t compt)
{
	auto & mesh_ { gradient.meshRef() };

	volumeField<vector> retField { mesh_, vector(0) };

	for (std::size_t i = 0; i < gradient.size(); ++i)
	{
		vector r;

		const vector & cellVector(mesh_.cells()[i].rC());

		std::vector<vector> neighbourGradientsOfSurfaces(
				mesh_.neighboursOfCells()[i].size());

		for (std::size_t cx = 0; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			const std::size_t cellIndex { mesh_.neighboursOfCells()[i][cx] };

			const vector deltaVector(
					mesh_.cells()[cellIndex].rC() - cellVector);

			const scalar deltaVectorMag { deltaVector.mag() };

			const vector deltaVectorNorm(deltaVector / deltaVectorMag);

			for (std::size_t j = 0; j < vector::vsize; ++j)
			{
				neighbourGradientsOfSurfaces[cx].r()[j] = (value()[cellIndex]
						- value()[i]) / deltaVectorMag * deltaVectorNorm()[j];
			}
		}

		if (mesh_.neighboursOfCells()[i].size()
				< mesh_.surfacesOfCells()[i].size())
		{
			for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
			{
				const std::size_t surfaceIndex { mesh_.surfacesOfCells()[i][j] };

				if (value.boundCond()[surfaceIndex].first
						!= boundaryConditionType::innerSurface)
				{
					vector deltaVector;
					if (value.boundCond()[surfaceIndex].first
							== boundaryConditionType::calculatedParallelBoundary)
						deltaVector =
								vector(
										(mesh_.surfaces()[surfaceIndex].rC()
												- cellVector)
												- parallel.cSdR().boundCond()[surfaceIndex].second);
					else
						deltaVector = vector(
								(mesh_.surfaces()[surfaceIndex].rC()
										- cellVector) * 2);

					const scalar deltaVectorMag { deltaVector.mag() };

					const vector deltaVectorNorm(deltaVector / deltaVectorMag);

					const scalar boundaryCellValue =
							bncCalc.boundaryConditionValueCell(value()[i],
									value.boundCond()[surfaceIndex],

									i, surfaceIndex, compt);

					neighbourGradientsOfSurfaces.emplace_back(
							vector(
									(boundaryCellValue - value()[i])
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									(boundaryCellValue - value()[i])
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									(boundaryCellValue - value()[i])
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm())));
				}
			}
		}
		else if (mesh_.neighboursOfCells()[i].size()
				> mesh_.surfacesOfCells()[i].size())
			throw exception(
					"Number of surfaces somehow lesser than number of neighbour cells.",
					errors::systemError);

		vector min(neighbourGradientsOfSurfaces[0]), max(
				neighbourGradientsOfSurfaces[0]);

		for (std::size_t cx = 1; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			for (std::size_t j = 0; j < vector::vsize; ++j)
				min.r()[j] = std::min(min()[j],
						neighbourGradientsOfSurfaces[cx]()[j]);

			for (std::size_t j = 0; j < vector::vsize; ++j)
				max.r()[j] = std::max(max()[j],
						neighbourGradientsOfSurfaces[cx]()[j]);
		}

		std::transform(min().begin(), min().end(), max().begin(), r.r().begin(),
				[](const auto num, const auto denom) 
				{	return num/(denom + stabilizator);});

		retField.r()[i] = limiterObjectP.calculate(r, gradient()[i]);
	}

	return retField;
}

schemi::volumeField<schemi::tensor> schemi::TVDLimiter(
		const volumeField<tensor> & gradient, const volumeField<vector> & value,
		const abstractLimiter & limiterObjectP,
		const boundaryConditionValue & bncCalc, const MPIHandler & parallel,
		const std::size_t compt)
{
	auto & mesh_ { gradient.meshRef() };

	volumeField<tensor> retField(mesh_, tensor(0));

	for (std::size_t i = 0; i < gradient.size(); ++i)
	{
		tensor r;

		const vector & cellVector(mesh_.cells()[i].rC());

		std::vector<tensor> neighbourGradientsOfSurfaces(
				mesh_.neighboursOfCells()[i].size());

		for (std::size_t cx = 0; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			const std::size_t cellIndex { mesh_.neighboursOfCells()[i][cx] };

			const vector deltaVector(
					mesh_.cells()[cellIndex].rC() - cellVector);

			const scalar deltaVectorMag { deltaVector.mag() };

			const vector deltaVectorNorm(deltaVector / deltaVectorMag);

			for (std::size_t j = 0; j < vector::vsize; ++j)
			{
				neighbourGradientsOfSurfaces[cx].r()[3 * j] = std::get<0>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 1] = std::get<1>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 2] = std::get<2>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
			}
		}

		if (mesh_.neighboursOfCells()[i].size()
				< mesh_.surfacesOfCells()[i].size())
		{
			for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
			{
				const std::size_t surfaceIndex = mesh_.surfacesOfCells()[i][j];

				if (value.boundCond()[surfaceIndex].first
						!= boundaryConditionType::innerSurface)
				{
					vector deltaVector;
					if (value.boundCond()[surfaceIndex].first
							== boundaryConditionType::calculatedParallelBoundary)
						deltaVector =
								vector(
										(mesh_.surfaces()[surfaceIndex].rC()
												- cellVector)
												- parallel.cSdR().boundCond()[surfaceIndex].second);
					else
						deltaVector = vector(
								(mesh_.surfaces()[surfaceIndex].rC()
										- cellVector) * 2);

					const scalar deltaVectorMag { deltaVector.mag() };

					const vector deltaVectorNorm(deltaVector / deltaVectorMag);

					const vector boundaryCellValue(
							bncCalc.boundaryConditionValueCell(value()[i],
									value.boundCond()[surfaceIndex], i,
									surfaceIndex, compt));

					neighbourGradientsOfSurfaces.emplace_back(
							tensor(
									std::get<0>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<1>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<2>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<0>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<1>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<2>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<0>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<1>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<2>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm())));
				}
			}
		}
		else if (mesh_.neighboursOfCells()[i].size()
				> mesh_.surfacesOfCells()[i].size())
			throw exception(
					"Number of surfaces somehow lesser than number of neighbour cells.",
					errors::systemError);

		tensor min(neighbourGradientsOfSurfaces[0]), max(
				neighbourGradientsOfSurfaces[0]);

		for (std::size_t cx = 1; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			for (std::size_t j = 0; j < tensor::vsize; ++j)
				min.r()[j] = std::min(min()[j],
						neighbourGradientsOfSurfaces[cx]()[j]);

			for (std::size_t j = 0; j < tensor::vsize; ++j)
				max.r()[j] = std::max(max()[j],
						neighbourGradientsOfSurfaces[cx]()[j]);
		}

		std::transform(min().begin(), min().end(), max().begin(), r.r().begin(),
				[](const auto num, const auto denom) 
				{	return num/(denom + stabilizator);});

		retField.r()[i] = limiterObjectP.calculate(r, gradient()[i]);
	}

	return retField;
}

schemi::volumeField<schemi::tensor3> schemi::TVDLimiter(
		const volumeField<tensor3> & gradient,
		const volumeField<tensor> & value,
		const abstractLimiter & limiterObjectP,
		const boundaryConditionValue & bncCalc, const MPIHandler & parallel,
		const std::size_t compt)
{
	auto & mesh_ { gradient.meshRef() };

	volumeField<tensor3> retField(mesh_, tensor3(0));

	for (std::size_t i = 0; i < gradient.size(); ++i)
	{
		tensor3 r;

		const vector & cellVector(mesh_.cells()[i].rC());

		std::vector<tensor3> neighbourGradientsOfSurfaces(
				mesh_.neighboursOfCells()[i].size());

		for (std::size_t cx = 0; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			const std::size_t cellIndex { mesh_.neighboursOfCells()[i][cx] };

			const vector deltaVector(
					mesh_.cells()[cellIndex].rC() - cellVector);

			const scalar deltaVectorMag { deltaVector.mag() };

			const vector deltaVectorNorm(deltaVector / deltaVectorMag);

			for (std::size_t j = 0; j < vector::vsize; ++j)
			{
				neighbourGradientsOfSurfaces[cx].r()[3 * j] = std::get<0>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 1] = std::get<1>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 2] = std::get<2>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 3] = std::get<3>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 4] = std::get<4>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 5] = std::get<5>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 6] = std::get<6>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 7] = std::get<7>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
				neighbourGradientsOfSurfaces[cx].r()[3 * j + 8] = std::get<8>(
						(value()[cellIndex] - value()[i])()) / deltaVectorMag
						* deltaVectorNorm()[j];
			}
		}

		if (mesh_.neighboursOfCells()[i].size()
				< mesh_.surfacesOfCells()[i].size())
		{
			for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
			{
				const std::size_t surfaceIndex = mesh_.surfacesOfCells()[i][j];

				if (value.boundCond()[surfaceIndex].first
						!= boundaryConditionType::innerSurface)
				{
					vector deltaVector;
					if (value.boundCond()[surfaceIndex].first
							== boundaryConditionType::calculatedParallelBoundary)
						deltaVector =
								vector(
										(mesh_.surfaces()[surfaceIndex].rC()
												- cellVector)
												- parallel.cSdR().boundCond()[surfaceIndex].second);
					else
						deltaVector = vector(
								(mesh_.surfaces()[surfaceIndex].rC()
										- cellVector) * 2);

					const scalar deltaVectorMag { deltaVector.mag() };

					const vector deltaVectorNorm(deltaVector / deltaVectorMag);

					const tensor boundaryCellValue(
							bncCalc.boundaryConditionValueCell(value()[i],
									value.boundCond()[surfaceIndex], i,
									surfaceIndex, compt));

					neighbourGradientsOfSurfaces.emplace_back(
							tensor3(
									std::get<0>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<1>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<2>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<3>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<4>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<5>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<6>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<7>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),
									std::get<8>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<0>(deltaVectorNorm()),

									std::get<0>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<1>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<2>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<3>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<4>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<5>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<6>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<7>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),
									std::get<8>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<1>(deltaVectorNorm()),

									std::get<0>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<1>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<2>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<3>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<4>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<5>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<6>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<7>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm()),
									std::get<8>(
											(boundaryCellValue - value()[i])())
											/ deltaVectorMag
											* std::get<2>(deltaVectorNorm())));
				}
			}
		}
		else if (mesh_.neighboursOfCells()[i].size()
				> mesh_.surfacesOfCells()[i].size())
			throw exception(
					"Number of surfaces somehow lesser than number of neighbour cells.",
					errors::systemError);

		tensor3 min(neighbourGradientsOfSurfaces[0]), max(
				neighbourGradientsOfSurfaces[0]);

		for (std::size_t cx = 1; cx < neighbourGradientsOfSurfaces.size(); ++cx)
		{
			for (std::size_t j = 0; j < tensor3::vsize; ++j)
				min.r()[j] = std::min(min()[j],
						neighbourGradientsOfSurfaces[cx]()[j]);

			for (std::size_t j = 0; j < tensor3::vsize; ++j)
				max.r()[j] = std::max(max()[j],
						neighbourGradientsOfSurfaces[cx]()[j]);
		}

		std::transform(min().begin(), min().end(), max().begin(), r.r().begin(),
				[](const auto num, const auto denom) 
				{	return num/(denom + stabilizator);});

		retField.r()[i] = limiterObjectP.calculate(r, gradient()[i]);
	}

	return retField;
}
