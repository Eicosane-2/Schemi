/*
 * instabilityParticlesHandler.hpp
 *
 *  Created on: 2025/08/25
 *      Author: Maxim Boldyrev
 */

#include <array>
#include <vector>

#include "boundaryConditionValue.hpp"
#include "BHRGoncharovTracerModel.hpp"
#include "MPIHandler.hpp"

#ifndef INSTABILITYPARTICLESHANDLER_HPP_
#define INSTABILITYPARTICLESHANDLER_HPP_

namespace schemi
{

class instabilityParticlesHandler
{
	constexpr static std::size_t maxRecursionDepthNearestCellSearch { 5 },
			maxRecursionDepthRhoSearch { 100 };

	constexpr static scalar purityParameter { 0.999 };

	enum positionType : std::size_t
	{
		rank = 0, cell, surface
	};

	const mesh & meshRef;
	const MPIHandler & parallelism;

	std::vector<std::size_t> boundarySurfaces;

	std::vector<std::array<std::size_t, 3>> particlePosition; // rank/cell/surface
	std::vector<std::array<scalar, 2>> cellSurfWeight;
	std::vector<BHRGoncharovTracerModel> particlesList;
	std::vector<initialisationStatus> particleStatus;
	std::size_t listSize { 0 };

	bool modelUsed;

	bool checkPosition(const vector & position) const noexcept;

	std::size_t findNearCell(const vector & position) noexcept;
	std::size_t findNearSurface(const vector & position,
			const std::size_t cellIndex) noexcept;
	std::array<scalar, 2> calculateWeights(const vector & position,
			const std::array<std::size_t, 3> & indexes);
	template<typename T>
	T calculateWeightedValue(const volumeField<T> & vCell,
			const surfaceField<T> & vSurf,
			const std::array<std::size_t, 3> & indexes,
			const std::array<schemi::scalar, 2> & cellSurfWeightParticle) const noexcept
	{
		const auto nearCellIndex = std::get<1>(indexes);
		const auto nearSurfIndex = std::get<2>(indexes);

		return vCell.cval()[nearCellIndex] * std::get<0>(cellSurfWeightParticle)
				+ vSurf.cval()[nearSurfIndex]
						* std::get<1>(cellSurfWeightParticle);
	}
	std::size_t searchOfNearestCell(const vector & newPosition,
			const std::pair<scalar, std::size_t> guessValue,
			std::vector<std::size_t> & checkedCells,
			std::size_t recursionLevel);
	std::array<std::size_t, 2> findNewNearPosition(const vector & newPosition,
			const std::size_t oldPositionCellIndex);

	void insertCaption(const std::size_t particleIndex) const;
	void clearLines(const std::size_t particleIndex,
			const std::pair<std::size_t, std::string> & readDataPoint) const;
	std::ifstream dataStream(const std::size_t particleIndex,
			const std::pair<std::size_t, std::string> & readDataPoint) const;

	void locateParticleNode(const std::size_t particleIndex,
			vector & positionVector,
			std::array<std::size_t, 1> & nodeParticleLocated);

	void searchDensity(const std::size_t sub,
			const std::array<std::size_t, 2> & mixSubs,
			const concentrationsPack<cubicCell> & concentrations,
			const std::vector<volumeField<scalar>> & densities,
			const vector & normale, const std::size_t cellIndex,
			scalar & rhoSub, vector & r, std::size_t recursionLevel,
			const boundaryConditionValue & boundVal,
			const std::valarray<scalar> & M) const;

	volumeField<scalar> formProfile(
			const std::array<std::size_t, 2> & substancesIndex,
			const vector & tracerPosition, const scalar radiusOfInfluence,
			const concentrationsPack<cubicCell> & concentrations,
			const boundaryConditionValue & bnc) const;

public:
	instabilityParticlesHandler(const mesh & meshIn, const MPIHandler & par,
			const volumeField<vector> & uCell,
			const surfaceField<vector> & uSurf,
			const std::pair<std::size_t, std::string> & readDataPoint);

	void timeIntegration(const volumeField<vector> & gradRho,
			const surfaceField<vector> & gradRhoSurf,
			const volumeField<vector> & uCell,
			const surfaceField<vector> & uSurf,
			const concentrationsPack<cubicCell> & concentrations,
			const std::vector<volumeField<scalar>> & densities,
			const boundaryConditionValue & boundVal,
			const std::valarray<scalar> & M, const scalar timestep,
			const volumeField<vector> & gradP, const volumeField<scalar> & divU,
			const volumeField<tensor> & gradU);

	void writeOutput(const std::string & fieldDataDirectoryName,
			const scalar Time) const;

	bool isInitialisationModelUsed() const noexcept
	{
		return modelUsed;
	}

	void checkTransitionToTurbulenceModel(const volumeField<scalar> & nuCell,
			const surfaceField<scalar> & nuSurface, volumeField<scalar> & k,
			volumeField<scalar> & epsilon, volumeField<vector> & a,
			volumeField<scalar> & b,
			const concentrationsPack<cubicCell> & concentrations,
			const boundaryConditionValue & boundVal, const scalar timestep);

	volumeField<scalar> generationArea(
			const concentrationsPack<cubicCell> & concentrations) const noexcept;
};
}

#endif /* INSTABILITYPARTICLESHANDLER_HPP_ */
