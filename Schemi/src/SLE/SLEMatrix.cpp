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

schemi::SLEMatrix::SLEMatrixStorage::SLEMatrixStorage(
		const mesh & meshRef) noexcept :
		centralDiagonale(0., meshRef.cellsSize()), freeTerm(0.,
				meshRef.cellsSize()), lowerTriangle(meshRef.cellsSize()), upperTriangle(
				meshRef.cellsSize())
{
}

void schemi::SLEMatrix::generateNabla(const volumeField<scalar> & vField,
		const volumeField<scalar> & rFieldOld,
		const volumeField<scalar> & rFieldNew,
		const surfaceField<vector> & additionalField, const scalar timestep,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	auto & mesh { vField.meshRef() };

	SLE.resize(1, SLEMatrixStorage(mesh));

	const scalar reversedTimestep { 1 / timestep };

	SLE[0].centralDiagonale = rFieldNew.ref() * reversedTimestep;

	SLE[0].freeTerm = rFieldOld.ref() * reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V();

				}
				else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceOwner()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V() * (-1);
				}
				else
					throw exception("Couldn't choose oIndex.",
							errorsEnum::systemError);

				if ((outerNormale & additionalField.ref()[surfaceIndex]) > 0)
				{
					SLE[0].centralDiagonale[i] +=
							(additionalField.ref()[surfaceIndex] & outerNormale);
				}
				else
				{
					if (oIndex < i)
						SLE[0].lowerTriangle[i].push_back(
								{ (additionalField.ref()[surfaceIndex]
										& outerNormale), oIndex });
					else
						SLE[0].upperTriangle[i].push_back(
								{ (additionalField.ref()[surfaceIndex]
										& outerNormale), oIndex });
				}
			}
				break;
			default:
			{
				const scalar boundaryValue { bncCalc.boundaryConditionValueCell(
						vField.ref()[i], vField.boundCond()[surfaceIndex], i,
						surfaceIndex, compt) };

				const vector outerNormale(
						mesh.surfaces()[surfaceIndex].N()
								* mesh.surfaces()[surfaceIndex].S()
								/ mesh.cells()[i].V());

				if ((outerNormale & additionalField.ref()[surfaceIndex]) > 0)
					SLE[0].centralDiagonale[i] +=
							(additionalField.ref()[surfaceIndex] & outerNormale);
				else
					SLE[0].freeTerm[i] -= (additionalField.ref()[surfaceIndex]
							& outerNormale) * boundaryValue;
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateLaplacian(const volumeField<scalar> & vField,
		const volumeField<scalar> & rFieldOld,
		const volumeField<scalar> & rFieldNew,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh { vField.meshRef() };

	SLE.resize(1, SLEMatrixStorage(mesh));

	const scalar reversedTimestep { 1 / timestep };

	SLE[0].centralDiagonale = rFieldNew.ref() * reversedTimestep;

	SLE[0].freeTerm = rFieldOld.ref() * reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V();
				}
				else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceOwner()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V() * (-1);
				}
				else
					throw exception("Couldn't choose oIndex.",
							errorsEnum::systemError);

				const vector deltaR(
						mesh.cells()[oIndex].rC() - mesh.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				if (oIndex < i)
					SLE[0].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
				else
					SLE[0].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
			}
				break;
			default:
			{
				const scalar boundaryValue { bncCalc.boundaryConditionValueCell(
						vField.ref()[i], vField.boundCond()[surfaceIndex], i,
						surfaceIndex, compt) };

				const vector deltaR(
						(mesh.surfaces()[surfaceIndex].rC()
								- mesh.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh.surfaces()[surfaceIndex].N()
								* mesh.surfaces()[surfaceIndex].S()
								/ mesh.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue;
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateLaplacian2TO(const volumeField<scalar> & vField,
		const volumeField<scalar> & rFieldOld,
		const volumeField<scalar> & rFieldNew,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh { vField.meshRef() };

	SLE.resize(1, SLEMatrixStorage(mesh));

	const scalar reversedTimestep { 1 / timestep };

	SLE[0].centralDiagonale = rFieldNew.ref() * reversedTimestep;

	SLE[0].freeTerm = rFieldOld.ref() * reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V();
				}
				else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceOwner()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V() * (-1);
				}
				else
					throw exception("Couldn't choose oIndex.",
							errorsEnum::systemError);

				const vector deltaR(
						mesh.cells()[oIndex].rC() - mesh.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i];

				if (oIndex < i)
					SLE[0].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
				else
					SLE[0].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex];
			}
				break;
			default:
			{
				const scalar boundaryValue { bncCalc.boundaryConditionValueCell(
						vField.ref()[i], vField.boundCond()[surfaceIndex], i,
						surfaceIndex, compt) };

				const vector deltaR(
						(mesh.surfaces()[surfaceIndex].rC()
								- mesh.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh.surfaces()[surfaceIndex].N()
								* mesh.surfaces()[surfaceIndex].S()
								/ mesh.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i];

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue;
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateLaplacian(const volumeField<vector> & vField,
		const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh { vField.meshRef() };

	SLE.resize(3, SLEMatrixStorage(mesh));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField.ref() * reversedTimestep;

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField.ref()[i] * vField.ref()[i].v()[j]
					* reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V();
				}
				else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceOwner()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V() * (-1);
				}
				else
					throw exception("Couldn't choose oIndex.",
							errorsEnum::systemError);

				const vector deltaR(
						mesh.cells()[oIndex].rC() - mesh.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				if (oIndex < i)
				{
					SLE[0].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[1].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[2].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
				}
				else
				{
					SLE[0].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[1].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[2].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
				}
			}
				break;
			default:
			{
				const vector boundaryValue(
						bncCalc.boundaryConditionValueCell(vField.ref()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh.surfaces()[surfaceIndex].rC()
								- mesh.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh.surfaces()[surfaceIndex].N()
								* mesh.surfaces()[surfaceIndex].S()
								/ mesh.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[0];
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[1];
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[2];
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateLaplacian2TO(const volumeField<vector> & vField,
		const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh { vField.meshRef() };

	SLE.resize(3, SLEMatrixStorage(mesh));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField.ref() * reversedTimestep;

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField.ref()[i] * vField.ref()[i].v()[j]
					* reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex = mesh.surfacesOfCells()[i][j];

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V();

				}
				else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceOwner()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V() * (-1);

				}
				else
					throw exception("Couldn't choose oIndex.",
							errorsEnum::systemError);

				const vector deltaR(
						mesh.cells()[oIndex].rC() - mesh.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[0];

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[1];

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[2];

				if (oIndex < i)
				{
					SLE[0].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[1].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[2].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
				}
				else
				{
					SLE[0].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[1].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[2].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
				}

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[0];
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[1];
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[2];
			}
				break;
			default:
			{
				const vector boundaryValue(
						bncCalc.boundaryConditionValueCell(vField.ref()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh.surfaces()[surfaceIndex].rC()
								- mesh.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh.surfaces()[surfaceIndex].N()
								* mesh.surfaces()[surfaceIndex].S()
								/ mesh.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[0];

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[1];

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[2];

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[0];
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[1];
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[2];
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateLaplacian(const volumeField<tensor> & vField,
		const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh { vField.meshRef() };

	SLE.resize(9, SLEMatrixStorage(mesh));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField.ref() * reversedTimestep;

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField.ref()[i] * vField.ref()[i].v()[j]
					* reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh.surfacesOfCells()[i][j] };

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V();

				}
				else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceOwner()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V() * (-1);

				}
				else
					throw exception("Couldn't choose oIndex.",
							errorsEnum::systemError);

				const vector deltaR(
						mesh.cells()[oIndex].rC() - mesh.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				if (oIndex < i)
				{
					SLE[0].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[1].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[2].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[3].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[4].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[5].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[6].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[7].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[8].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
				}
				else
				{
					SLE[0].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[1].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[2].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[3].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[4].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[5].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[6].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[7].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
					SLE[8].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm),
									oIndex });
				}
			}
				break;
			default:
			{
				const tensor boundaryValue(
						bncCalc.boundaryConditionValueCell(vField.ref()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh.surfaces()[surfaceIndex].rC()
								- mesh.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh.surfaces()[surfaceIndex].N()
								* mesh.surfaces()[surfaceIndex].S()
								/ mesh.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[0];
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[1];
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[2];
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[3];
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[4];
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[5];
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[6];
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[7];
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[8];
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateLaplacian2TO(const volumeField<tensor> & vField,
		const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const scalar timestep, const boundaryConditionValue & bncCalc,
		const std::size_t compt)
{
	auto & mesh { vField.meshRef() };

	SLE.resize(9, SLEMatrixStorage(mesh));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField.ref() * reversedTimestep;

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField.ref()[i] * vField.ref()[i].v()[j]
					* reversedTimestep;

	for (std::size_t i = 0; i < vField.size(); ++i)
		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex = mesh.surfacesOfCells()[i][j];

			switch (vField.boundCond()[surfaceIndex].first)
			{
			case boundaryConditionType::innerSurface:
			{
				std::size_t oIndex;

				vector outerNormale;

				if (mesh.surfaceOwner()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceNeighbour()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V();

				}
				else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				{
					oIndex = mesh.surfaceOwner()[surfaceIndex];

					outerNormale = mesh.surfaces()[surfaceIndex].N()
							* mesh.surfaces()[surfaceIndex].S()
							/ mesh.cells()[i].V() * (-1);

				}
				else
					throw exception("Couldn't choose oIndex.",
							errorsEnum::systemError);

				const vector deltaR(
						mesh.cells()[oIndex].rC() - mesh.cells()[i].rC());

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[0];

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[1];

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[2];

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[3];

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[4];

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[5];

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[6];

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[7];

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[8];

				if (oIndex < i)
				{
					SLE[0].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[1].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[2].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[3].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[4].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[5].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[6].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[7].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[8].lowerTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
				}
				else
				{
					SLE[0].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[1].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[2].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[3].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[4].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[5].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[6].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[7].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
					SLE[8].upperTriangle[i].push_back(
							{ -effectiveDiffusionCoefficient.ref()[surfaceIndex]
									/ deltaRMag * (outerNormale & deltaRNorm)
									* 0.5, oIndex });
				}

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[0];
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[1];
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[2];
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[3];
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[4];
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[5];
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[6];
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[7];
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[oIndex].v()[8];
			}
				break;
			default:
			{
				const tensor boundaryValue(
						bncCalc.boundaryConditionValueCell(vField.ref()[i],
								vField.boundCond()[surfaceIndex], i,
								surfaceIndex, compt));

				const vector deltaR(
						(mesh.surfaces()[surfaceIndex].rC()
								- mesh.cells()[i].rC()) * 2);

				const scalar deltaRMag { deltaR.mag() };

				const vector deltaRNorm(deltaR / deltaRMag);

				const vector outerNormale(
						mesh.surfaces()[surfaceIndex].N()
								* mesh.surfaces()[surfaceIndex].S()
								/ mesh.cells()[i].V());

				SLE[0].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[0];

				SLE[1].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm);
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[1];

				SLE[2].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[2];

				SLE[3].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[3];

				SLE[4].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[4];

				SLE[5].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[5];

				SLE[6].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[6];

				SLE[7].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[7];

				SLE[8].centralDiagonale[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5;
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm) * 0.5
								* vField.ref()[i].v()[8];

				SLE[0].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[0];
				SLE[1].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[1];
				SLE[2].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[2];
				SLE[3].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[3];
				SLE[4].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[4];
				SLE[5].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[5];
				SLE[6].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[6];
				SLE[7].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[7];
				SLE[8].freeTerm[i] +=
						effectiveDiffusionCoefficient.ref()[surfaceIndex]
								/ deltaRMag * (outerNormale & deltaRNorm)
								* boundaryValue.v()[8];
			}
				break;
			}
		}
}

void schemi::SLEMatrix::generateExplicitLaplacian(
		const volumeField<scalar> & vField,
		const volumeField<scalar> & rFieldOld,
		const volumeField<scalar> & rFieldNew,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const surfaceField<scalar> & surfaceVField, const scalar timestep,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	SLE.resize(1, SLEMatrixStorage(vField.meshRef()));

	const scalar reversedTimestep { 1 / timestep };

	SLE[0].centralDiagonale = rFieldNew.ref();

	SLE[0].freeTerm = rFieldOld.ref() * reversedTimestep;

	auto gradInCell = grad(surfaceVField);

	auto flowThroughSurface = gradientLinearInterpolate(gradInCell, vField,
			bncCalc, compt);

	for (std::size_t i = 0; i < flowThroughSurface.size(); ++i)
		flowThroughSurface.ref_r()[i] = flowThroughSurface.ref()[i]
				* effectiveDiffusionCoefficient.ref()[i] * (-1);

	auto divFlow = divergence(flowThroughSurface);

	SLE[0].freeTerm -= divFlow.ref();
}

void schemi::SLEMatrix::generateExplicitLaplacian(
		const volumeField<vector> & vField, const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const surfaceField<vector> & surfaceVField, const scalar timestep,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	SLE.resize(3, SLEMatrixStorage(vField.meshRef()));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField.ref();

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField.ref()[i] * vField.ref()[i].v()[j]
					* reversedTimestep;

	auto gradInCell = grad(surfaceVField);

	auto flowThroughSurface = gradientLinearInterpolate(gradInCell, vField,
			bncCalc, compt);

	for (std::size_t i = 0; i < flowThroughSurface.size(); ++i)
		flowThroughSurface.ref_r()[i] = flowThroughSurface.ref()[i]
				* effectiveDiffusionCoefficient.ref()[i] * (-1);

	auto divFlow = divergence(flowThroughSurface);

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] -= divFlow.ref()[i].v()[j];
}

void schemi::SLEMatrix::generateExplicitLaplacian(
		const volumeField<tensor> & vField, const volumeField<scalar> & rField,
		const surfaceField<scalar> & effectiveDiffusionCoefficient,
		const surfaceField<tensor> & surfaceVField, const scalar timestep,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	SLE.resize(9, SLEMatrixStorage(vField.meshRef()));

	const scalar reversedTimestep { 1 / timestep };

	for (std::size_t j = 0; j < SLE.size(); ++j)
		SLE[j].centralDiagonale = rField.ref();

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] = rField.ref()[i] * vField.ref()[i].v()[j]
					* reversedTimestep;

	auto gradInCell = grad(surfaceVField);

	auto flowThroughSurface = gradientLinearInterpolate(gradInCell, vField,
			bncCalc, compt);

	for (std::size_t i = 0; i < flowThroughSurface.size(); ++i)
		flowThroughSurface.ref_r()[i] = flowThroughSurface.ref()[i]
				* effectiveDiffusionCoefficient.ref()[i] * (-1);

	auto divFlow = divergence(flowThroughSurface);

	for (std::size_t j = 0; j < SLE.size(); ++j)
		for (std::size_t i = 0; i < vField.size(); ++i)
			SLE[j].freeTerm[i] -= divFlow.ref()[i].v()[j];
}

void schemi::SLEMatrix::distributeMinusSourceTerm(
		const volumeField<scalar> & source,
		const volumeField<scalar> & basicField) noexcept
{
	for (std::size_t i = 0; i < source.size(); ++i)
		if (source.ref()[i] <= 0.)
			SLE[0].freeTerm[i] -= source.ref()[i];
		else
			SLE[0].centralDiagonale[i] += source.ref()[i]
					/ (basicField.ref()[i] + stabilizator);
}

void schemi::SLEMatrix::distributePlusSourceTerm(
		const volumeField<scalar> & source,
		const volumeField<scalar> & basicField) noexcept
{
	for (std::size_t i = 0; i < source.size(); ++i)
		if (source.ref()[i] >= 0.)
			SLE[0].freeTerm[i] += source.ref()[i];
		else
			SLE[0].centralDiagonale[i] -= source.ref()[i]
					/ (basicField.ref()[i] + stabilizator);
}
