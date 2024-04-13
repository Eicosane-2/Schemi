/*
 * SLEMatrix.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "SLEMatrix.hpp"

#include "vectorVectorDotProduct.hpp"
#include "divergence.hpp"
#include "gradient.hpp"

schemi::SLEMatrix::SLEMatrix(const std::string & stringIn) noexcept :
		name(stringIn)
{
}

schemi::SLEMatrix::SLEMatrixStorage::SLEMatrixStorage() :
		centralDiagonale(0), freeTerm(0), lowerTriangle(0), upperTriangle(0)
{
}

schemi::SLEMatrix::SLEMatrixStorage::SLEMatrixStorage(
		const mesh & meshRef) noexcept :
		centralDiagonale(0., meshRef.cellsSize()), freeTerm(0.,
				meshRef.cellsSize()), lowerTriangle(meshRef.cellsSize()), upperTriangle(
				meshRef.cellsSize())
{
}

void schemi::SLEMatrix::SLEMatrixStorage::transpose() noexcept
{
	std::vector<std::vector<std::pair<scalar, std::size_t>>> lowerTriangleNew(
			lowerTriangle.size()), upperTriangleNew(upperTriangle.size());

	for (std::size_t i = 0; i < centralDiagonale.size(); ++i)
	{
		for (std::size_t j = 0; j < lowerTriangle[i].size(); ++j)
		{
			const std::size_t jAbsOld = lowerTriangle[i][j].second;
			const auto Avalue = lowerTriangle[i][j].first;

			const std::size_t iNew = jAbsOld;
			const std::size_t jNew = i;

			upperTriangleNew[iNew].emplace_back(Avalue, jNew);
		}

		for (std::size_t j = 0; j < upperTriangle[i].size(); ++j)
		{
			const std::size_t jAbsOld = upperTriangle[i][j].second;
			const auto Avalue = upperTriangle[i][j].first;

			const std::size_t iNew = jAbsOld;
			const std::size_t jNew = i;

			lowerTriangleNew[iNew].emplace_back(Avalue, jNew);
		}
	}

	lowerTriangle = lowerTriangleNew;
	upperTriangle = upperTriangleNew;
}

void schemi::SLEMatrix::generateLaplacianSurfaceBoundary(
		const volumeField<scalar> & vField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	auto & mesh_ { vField.meshRef() };

	SLE.resize(1, SLEMatrixStorage(mesh_));

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh_.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V();
				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceOwner()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V() * (-1);
				}
				else
					throw exception("Couldn't choose oIndex.",
							errors::systemError);

				const vector deltaR(
						mesh_.cells()[oIndex].rC() - mesh_.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				if (oIndex < i)
					SLE[0].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
				else
					SLE[0].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
			}
				break;
			case boundaryConditionType::calculatedParallelBoundary:
			{
				const scalar boundaryValue { bncCalc.boundaryConditionValueCell(
						vField()[i], vField.boundCond()[surfaceIndex], i,
						surfaceIndex, compt) };

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC())
								- bncCalc.parallelism.cSdR().boundCond()[surfaceIndex].second);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue;
			}
				break;
			default:
			{
				const scalar boundaryValue {
						bncCalc.boundaryConditionValueSurface(vField()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt) };

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC()));

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue;
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateDTimeNabla(const volumeField<scalar> & vField,
		const volumeField<scalar> & rFieldOld,
		const volumeField<scalar> & rFieldNew,
		const surfaceField<vector> & additionalField, const scalar timestep,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	auto & mesh_ { vField.meshRef() };

	SLE.resize(1, SLEMatrixStorage(mesh_));

	const scalar reversedTimestep { 1 / timestep };

	SLE[0].centralDiagonale = rFieldNew() * reversedTimestep;

	SLE[0].freeTerm = rFieldOld() * reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh_.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V();

				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceOwner()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V() * (-1);
				}
				else
					throw exception("Couldn't choose oIndex.",
							errors::systemError);

				if ((outerNormale & additionalField()[surfaceIndex]) > 0)
				{
					SLE[0].centralDiagonale[i] +=
							(additionalField()[surfaceIndex] & outerNormale);
				}
				else
				{
					if (oIndex < i)
						SLE[0].lowerTriangle[i].emplace_back(
								(additionalField()[surfaceIndex] & outerNormale),
								oIndex);
					else
						SLE[0].upperTriangle[i].emplace_back(
								(additionalField()[surfaceIndex] & outerNormale),
								oIndex);
				}
			}
				break;
			default:
			{
				const scalar boundaryValue { bncCalc.boundaryConditionValueCell(
						vField()[i], vField.boundCond()[surfaceIndex], i,
						surfaceIndex, compt) };

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				if ((outerNormale & additionalField()[surfaceIndex]) > 0)
					SLE[0].centralDiagonale[i] +=
							(additionalField()[surfaceIndex] & outerNormale);
				else
					SLE[0].freeTerm[i] -= (additionalField()[surfaceIndex]
							& outerNormale) * boundaryValue;
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateDTimeLaplacian(
		const volumeField<scalar> & vField,
		const volumeField<scalar> & rFieldOld,
		const volumeField<scalar> & rFieldNew,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh_ { vField.meshRef() };

	SLE.resize(1, SLEMatrixStorage(mesh_));

	const scalar reversedTimestep { 1 / timestep };

	SLE[0].centralDiagonale = rFieldNew() * reversedTimestep;

	SLE[0].freeTerm = rFieldOld() * reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh_.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V();
				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceOwner()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V() * (-1);
				}
				else
					throw exception("Couldn't choose oIndex.",
							errors::systemError);

				const vector deltaR(
						mesh_.cells()[oIndex].rC() - mesh_.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				if (oIndex < i)
					SLE[0].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
				else
					SLE[0].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
			}
				break;
			case boundaryConditionType::calculatedParallelBoundary:
			{
				const scalar boundaryValue { bncCalc.boundaryConditionValueCell(
						vField()[i], vField.boundCond()[surfaceIndex], i,
						surfaceIndex, compt) };

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC())
								- bncCalc.parallelism.cSdR().boundCond()[surfaceIndex].second);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue;
			}
				break;
			default:
			{
				const scalar boundaryValue { bncCalc.boundaryConditionValueCell(
						vField()[i], vField.boundCond()[surfaceIndex], i,
						surfaceIndex, compt) };

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue;
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateDTimeLaplacian2TO(
		const volumeField<scalar> & vField,
		const volumeField<scalar> & rFieldOld,
		const volumeField<scalar> & rFieldNew,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh_ { vField.meshRef() };

	SLE.resize(1, SLEMatrixStorage(mesh_));

	const scalar reversedTimestep { 1 / timestep };

	SLE[0].centralDiagonale = rFieldNew() * reversedTimestep;

	SLE[0].freeTerm = rFieldOld() * reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh_.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V();
				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceOwner()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V() * (-1);
				}
				else
					throw exception("Couldn't choose oIndex.",
							errors::systemError);

				const vector deltaR(
						mesh_.cells()[oIndex].rC() - mesh_.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField()[i];

				if (oIndex < i)
					SLE[0].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
				else
					SLE[0].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField()[oIndex];
			}
				break;
			case boundaryConditionType::calculatedParallelBoundary:
			{
				const scalar boundaryValue { bncCalc.boundaryConditionValueCell(
						vField()[i], vField.boundCond()[surfaceIndex], i,
						surfaceIndex, compt) };

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC())
								- bncCalc.parallelism.cSdR().boundCond()[surfaceIndex].second);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField()[i];

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue;
			}
				break;
			default:
			{
				const scalar boundaryValue { bncCalc.boundaryConditionValueCell(
						vField()[i], vField.boundCond()[surfaceIndex], i,
						surfaceIndex, compt) };

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField()[i];

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue;
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateDTimeLaplacian(
		const volumeField<vector> & vField, const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh_ { vField.meshRef() };

	SLE.resize(3, SLEMatrixStorage(mesh_));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField() * reversedTimestep;

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField()[i] * vField()[i]()[j]
					* reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh_.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V();
				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceOwner()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V() * (-1);
				}
				else
					throw exception("Couldn't choose oIndex.",
							errors::systemError);

				const vector deltaR(
						mesh_.cells()[oIndex].rC() - mesh_.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				if (oIndex < i)
				{
					SLE[0].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[1].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[2].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
				}
				else
				{
					SLE[0].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[1].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[2].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
				}
			}
				break;
			case boundaryConditionType::calculatedParallelBoundary:
			{
				const vector boundaryValue(
						bncCalc.boundaryConditionValueCell(vField()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC())
								- bncCalc.parallelism.cSdR().boundCond()[surfaceIndex].second);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<0>(boundaryValue());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<1>(boundaryValue());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<2>(boundaryValue());
			}
				break;
			default:
			{
				const vector boundaryValue(
						bncCalc.boundaryConditionValueCell(vField()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<0>(boundaryValue());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<1>(boundaryValue());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<2>(boundaryValue());
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateDTimeLaplacian2TO(
		const volumeField<vector> & vField, const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh_ { vField.meshRef() };

	SLE.resize(3, SLEMatrixStorage(mesh_));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField() * reversedTimestep;

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField()[i] * vField()[i]()[j]
					* reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex = mesh_.surfacesOfCells()[i][j];

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V();

				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceOwner()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V() * (-1);

				}
				else
					throw exception("Couldn't choose oIndex.",
							errors::systemError);

				const vector deltaR(
						mesh_.cells()[oIndex].rC() - mesh_.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<0>(vField()[i]());

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<1>(vField()[i]());

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<2>(vField()[i]());

				if (oIndex < i)
				{
					SLE[0].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[1].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[2].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
				}
				else
				{
					SLE[0].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[1].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[2].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
				}

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<0>(vField()[oIndex]());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<1>(vField()[oIndex]());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<2>(vField()[oIndex]());
			}
				break;
			case boundaryConditionType::calculatedParallelBoundary:
			{
				const vector boundaryValue(
						bncCalc.boundaryConditionValueCell(vField()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC())
								- bncCalc.parallelism.cSdR().boundCond()[surfaceIndex].second);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<0>(vField()[i]());

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<1>(vField()[i]());

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<2>(vField()[i]());

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<0>(boundaryValue());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<1>(boundaryValue());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<2>(boundaryValue());
			}
				break;
			default:
			{
				const vector boundaryValue(
						bncCalc.boundaryConditionValueCell(vField()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<0>(vField()[i]());

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<1>(vField()[i]());

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<2>(vField()[i]());

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<0>(boundaryValue());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<1>(boundaryValue());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<2>(boundaryValue());
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateDTimeLaplacian(
		const volumeField<tensor> & vField, const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh_ { vField.meshRef() };

	SLE.resize(9, SLEMatrixStorage(mesh_));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField() * reversedTimestep;

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField()[i] * vField()[i]()[j]
					* reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh_.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V();

				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceOwner()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V() * (-1);

				}
				else
					throw exception("Couldn't choose oIndex.",
							errors::systemError);

				const vector deltaR(
						mesh_.cells()[oIndex].rC() - mesh_.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				if (oIndex < i)
				{
					SLE[0].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[1].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[2].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[3].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[4].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[5].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[6].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[7].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[8].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
				}
				else
				{
					SLE[0].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[1].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[2].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[3].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[4].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[5].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[6].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[7].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
					SLE[8].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
							oIndex);
				}
			}
				break;
			case boundaryConditionType::calculatedParallelBoundary:
			{
				const tensor boundaryValue(
						bncCalc.boundaryConditionValueCell(vField()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC())
								- bncCalc.parallelism.cSdR().boundCond()[surfaceIndex].second);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<0>(boundaryValue());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<1>(boundaryValue());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<2>(boundaryValue());
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<3>(boundaryValue());
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<4>(boundaryValue());
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<5>(boundaryValue());
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<6>(boundaryValue());
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<7>(boundaryValue());
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<8>(boundaryValue());
			}
				break;
			default:
			{
				const tensor boundaryValue(
						bncCalc.boundaryConditionValueCell(vField()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<0>(boundaryValue());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<1>(boundaryValue());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<2>(boundaryValue());
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<3>(boundaryValue());
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<4>(boundaryValue());
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<5>(boundaryValue());
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<6>(boundaryValue());
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<7>(boundaryValue());
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<8>(boundaryValue());
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateDTimeLaplacian2TO(
		const volumeField<tensor> & vField, const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh_ { vField.meshRef() };

	SLE.resize(9, SLEMatrixStorage(mesh_));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField() * reversedTimestep;

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField()[i] * vField()[i]()[j]
					* reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex = mesh_.surfacesOfCells()[i][j];

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V();

				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh_.surfaceOwner()[surfaceIndex];

					outerNormale = mesh_.surfaces()[surfaceIndex].N()
							* mesh_.surfaces()[surfaceIndex].S()
							/ mesh_.cells()[i].V() * (-1);

				}
				else
					throw exception("Couldn't choose oIndex.",
							errors::systemError);

				const vector deltaR(
						mesh_.cells()[oIndex].rC() - mesh_.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<0>(vField()[i]());

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<1>(vField()[i]());

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<2>(vField()[i]());

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<3>(vField()[i]());

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<4>(vField()[i]());

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<5>(vField()[i]());

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<6>(vField()[i]());

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<7>(vField()[i]());

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<8>(vField()[i]());

				if (oIndex < i)
				{
					SLE[0].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[1].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[2].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[3].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[4].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[5].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[6].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[7].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[8].lowerTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
				}
				else
				{
					SLE[0].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[1].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[2].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[3].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[4].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[5].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[6].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[7].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
					SLE[8].upperTriangle[i].emplace_back(
							-effectiveDiffusionCoefficient()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex);
				}

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<0>(vField()[oIndex]());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<1>(vField()[oIndex]());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<2>(vField()[oIndex]());
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<3>(vField()[oIndex]());
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<4>(vField()[oIndex]());
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<5>(vField()[oIndex]());
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<6>(vField()[oIndex]());
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<7>(vField()[oIndex]());
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<8>(vField()[oIndex]());
			}
				break;
			case boundaryConditionType::calculatedParallelBoundary:
			{
				const tensor boundaryValue(
						bncCalc.boundaryConditionValueCell(vField()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC())
								- bncCalc.parallelism.cSdR().boundCond()[surfaceIndex].second);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<0>(vField()[i]());

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<1>(vField()[i]());

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<2>(vField()[i]());

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<3>(vField()[i]());

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<4>(vField()[i]());

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<5>(vField()[i]());

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<6>(vField()[i]());

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<7>(vField()[i]());

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<8>(vField()[i]());

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<0>(boundaryValue());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<1>(boundaryValue());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<2>(boundaryValue());
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<3>(boundaryValue());
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<4>(boundaryValue());
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<5>(boundaryValue());
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<6>(boundaryValue());
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<7>(boundaryValue());
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<8>(boundaryValue());
			}
				break;
			default:
			{
				const tensor boundaryValue(
						bncCalc.boundaryConditionValueCell(vField()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh_.surfaces()[surfaceIndex].rC()
								- mesh_.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh_.surfaces()[surfaceIndex].N()
								* mesh_.surfaces()[surfaceIndex].S()
								/ mesh_.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<0>(vField()[i]());

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<1>(vField()[i]());

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<2>(vField()[i]());

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<3>(vField()[i]());

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<4>(vField()[i]());

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<5>(vField()[i]());

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<6>(vField()[i]());

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<7>(vField()[i]());

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* std::get<8>(vField()[i]());

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<0>(boundaryValue());
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<1>(boundaryValue());
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<2>(boundaryValue());
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<3>(boundaryValue());
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<4>(boundaryValue());
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<5>(boundaryValue());
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<6>(boundaryValue());
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<7>(boundaryValue());
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* std::get<8>(boundaryValue());
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateDTimeExplicitLaplacian(
		const volumeField<scalar> & vField,
		const volumeField<scalar> & rFieldOld,
		const volumeField<scalar> & rFieldNew,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	SLE.resize(1, SLEMatrixStorage(vField.meshRef()));

	SLE[0].centralDiagonale = rFieldNew();

	SLE[0].explOldTime.resize(SLE[0].freeTerm.size(), 0.0);

	SLE[0].explOldTime = rFieldOld();

	auto flowThroughSurface = surfGrad(vField, bncCalc, compt);

	for (std::size_t i = 0; i < flowThroughSurface.size(); ++i)
		flowThroughSurface.r()[i] = flowThroughSurface()[i]
				* effectiveDiffusionCoefficient()[i] * (-1);

	auto divFlow = divergence(flowThroughSurface);

	SLE[0].freeTerm -= divFlow();
}

void schemi::SLEMatrix::generateDTimeExplicitLaplacian(
		const volumeField<vector> & vField, const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	SLE.resize(3, SLEMatrixStorage(vField.meshRef()));

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField();

	for (std::size_t j = 0; j < SLE.size(); ++j)
	{
		SLE[j].explOldTime.resize(SLE[j].freeTerm.size(), 0.0);
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].explOldTime[i] = rField()[i] * vField()[i]()[j];
	}

	auto flowThroughSurface = surfGrad(vField, bncCalc, compt);

	for (std::size_t i = 0; i < flowThroughSurface.size(); ++i)
		flowThroughSurface.r()[i] = flowThroughSurface()[i]
				* effectiveDiffusionCoefficient()[i] * (-1);

	auto divFlow = divergence(flowThroughSurface);

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] -= divFlow()[i]()[j];
}

void schemi::SLEMatrix::generateDTimeExplicitLaplacian(
		const volumeField<tensor> & vField, const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	SLE.resize(9, SLEMatrixStorage(vField.meshRef()));

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField();

	for (std::size_t j = 0; j < SLE.size(); ++j)
	{
		SLE[j].explOldTime.resize(SLE[j].freeTerm.size(), 0.0);
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].explOldTime[i] = rField()[i] * vField()[i]()[j];
	}

	auto flowThroughSurface = surfGrad(vField, bncCalc, compt);

	for (std::size_t i = 0; i < flowThroughSurface.size(); ++i)
		flowThroughSurface.r()[i] = flowThroughSurface()[i]
				* effectiveDiffusionCoefficient()[i] * (-1);

	auto divFlow = divergence(flowThroughSurface);

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] -= divFlow()[i]()[j];
}

void schemi::SLEMatrix::freeSourceTerm(
		const volumeField<scalar> & source) noexcept
{
	SLE[0].freeTerm += source();
}

void schemi::SLEMatrix::freeSourceTerm(
		const volumeField<vector> & source) noexcept
{
	for (std::size_t j = 0; j < vector::vsize; j++)
		for (std::size_t i = 0; i < source.size(); ++i)
			SLE[j].freeTerm[i] += source()[i]()[j];
}

void schemi::SLEMatrix::freeSourceTerm(
		const volumeField<tensor> & source) noexcept
{
	for (std::size_t j = 0; j < tensor::vsize; j++)
		for (std::size_t i = 0; i < source.size(); ++i)
			SLE[j].freeTerm[i] += source()[i]()[j];
}

void schemi::SLEMatrix::distributeSourceTerm(const volumeField<scalar> & source,
		const volumeField<scalar> & basicField) noexcept
{
	for (std::size_t i = 0; i < source.size(); ++i)
		if (source()[i] >= 0.)
			SLE[0].freeTerm[i] += source()[i];
		else
			SLE[0].centralDiagonale[i] -= source()[i]
					/ (basicField()[i] + stabilizator);
}

void schemi::SLEMatrix::distributeSourceTerm(const volumeField<vector> & source,
		const volumeField<vector> & basicField) noexcept
{
	for (std::size_t j = 0; j < vector::vsize; j++)
		for (std::size_t i = 0; i < source.size(); ++i)
			if (source()[i]()[j] >= 0.)
				SLE[j].freeTerm[i] += source()[i]()[j];
			else
				SLE[j].centralDiagonale[i] -= source()[i]()[j]
						/ (basicField()[i]()[j] + stabilizator);
}

void schemi::SLEMatrix::distributeSourceTerm(const volumeField<tensor> & source,
		const volumeField<tensor> & basicField) noexcept
{
	for (std::size_t j = 0; j < tensor::vsize; j++)
		for (std::size_t i = 0; i < source.size(); ++i)
			if (source()[i]()[j] >= 0.)
				SLE[j].freeTerm[i] += source()[i]()[j];
			else
				SLE[j].centralDiagonale[i] -= source()[i]()[j]
						/ (basicField()[i]()[j] + stabilizator);
}

void schemi::SLEMatrix::diagonaleSourceTerm(
		const volumeField<scalar> & source) noexcept
{
	SLE[0].centralDiagonale -= source();
}

void schemi::SLEMatrix::diagonaleSourceTerm(
		const volumeField<vector> & source) noexcept
{
	for (std::size_t j = 0; j < vector::vsize; j++)
		for (std::size_t i = 0; i < source.size(); ++i)
			SLE[j].centralDiagonale[i] -= source()[i]()[j];
}

void schemi::SLEMatrix::diagonaleSourceTerm(
		const volumeField<tensor> & source) noexcept
{
	for (std::size_t j = 0; j < tensor::vsize; j++)
		for (std::size_t i = 0; i < source.size(); ++i)
			SLE[j].centralDiagonale[i] -= source()[i]()[j];
}
