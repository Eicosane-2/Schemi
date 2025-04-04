/*
 * MPIHandler.hpp
 *
 *  Created on: 2020/09/07
 *      Author: Maxim Boldyrev
 *
 *      Functions for MPI data exchange.
 */

#ifndef MPIHANDLER_HPP_
#define MPIHANDLER_HPP_

#include <cstddef>
#include <memory>

#include "scalar.hpp"
#include "tensor.hpp"
#include "tensor3.hpp"
#include "vector.hpp"
#include "bunchOfFields.hpp"
#include "surfaceField.hpp"
#include "volumeField.hpp"

#if defined(MPI_ENABLE) && defined(MPI_DEBUG)
#define MPI_VERSION 111
#define MPI_COMM_WORLD 222
#define MPI_UNSIGNED 202
#define MPI_UNSIGNED_LONG 303
#define MPI_UNSIGNED_LONG_LONG 404
#define MPI_DOUBLE 505
#define MPI_STATUSES_IGNORE 123

int MPI_Init(int*, char***);
int MPI_Comm_rank(int, int*);
int MPI_Comm_size(int, int*);
int MPI_Barrier(int);
int MPI_Gather(const std::size_t*, int, int, std::size_t*, int, int, int, int);
int MPI_Gather(const double*, int, int, double*, int, int, int, int);
int MPI_Send(const double*, int, int, int, int, int);
int MPI_Recv(double*, int, int, int, int, int, int);
int MPI_Bcast(double*, int, int, int, int);
int MPI_Bcast(std::size_t*, int, int, int, int);
int MPI_Finalize();
#elif defined(MPI_ENABLE)
#include <mpi.h>
#endif

namespace schemi
{
#if SIZE_MAX == UINT_MAX
#define schemi_MPI_SIZE MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define schemi_MPI_SIZE MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define schemi_MPI_SIZE MPI_UNSIGNED_LONG_LONG
#else
#error schemi_MPI_SIZE undefined.
#endif

class MPIHandler
{
public:
	typedef double mpi_scalar;

private:
	constexpr static std::size_t tailSurfaceStart = 0;
	std::size_t tailSurfaceEnd = 0,

	pointSurfaceStart = 0, pointSurfaceEnd = 0,

	bottomSurfaceStart = 0, bottomSurfaceEnd = 0,

	rightSurfaceStart = 0, rightSurfaceEnd = 0,

	leftSurfaceStart = 0, leftSurfaceEnd = 0,

	topSurfaceStart = 0, topSurfaceEnd = 0,

	gathBufSize = 0, localCellSize[1] { 0 };

	std::unique_ptr<mpi_scalar> tailData_to_send { nullptr },
			pointData_to_send { nullptr }, bottomData_to_send { nullptr },
			rightData_to_send { nullptr }, leftData_to_send { nullptr },
			topData_to_send { nullptr }, tailData_to_rec { nullptr },
			pointData_to_rec { nullptr }, bottomData_to_rec { nullptr },
			rightData_to_rec { nullptr }, leftData_to_rec { nullptr },
			topData_to_rec { nullptr },

			gathBuf { nullptr }, cellSendBuf { nullptr };

	std::array<std::size_t, 3> rankDir { 0, 0, 0 }, rankDirMax { 0, 0, 0 };

	bool arraysInit = false;

	enum class MPITags
	{
		tail, point, bottom, right, left, top
	};

	struct
	{
		int & tail { a[0] }, & point { a[1] }, & bottom { a[2] },
				& right { a[3] }, & left { a[4] }, & top { a[5] };
	private:
		std::array<int, 6> a { -1, -1, -1, -1, -1, -1 };
	} envNodes;

public:
	const std::size_t mpi_rank, mpi_size, root = 0;

	std::unique_ptr<volumeField<scalar>> parallCellVolume { nullptr };
	std::unique_ptr<surfaceField<vector>> ownerSurfaceDeltaR { nullptr };

	std::size_t totCellNum() const noexcept;

	const volumeField<scalar>& Vol() const;

	const surfaceField<vector>& cSdR() const;

	MPIHandler(std::size_t mpi_rank_in, std::size_t mpi_size_in);

	MPIHandler(const MPIHandler&) = delete;
	auto& operator=(const MPIHandler&) = delete;

	std::pair<vector, vector> correctParallelepipedVector(
			const vector & parallelepiped) const noexcept;

	void initialiseBuffersSize([[maybe_unused]] const mesh & mesh_);

	void initialiseParalleMeshData(const mesh & mesh_);

	void correctBoundaryValues(
			[[maybe_unused]] volumeField<scalar> & corField) const noexcept;

	void correctBoundaryValues(
			[[maybe_unused]] volumeField<vector> & corField) const noexcept;

	void correctBoundaryValues(
			[[maybe_unused]] volumeField<tensor> & corField) const noexcept;

	void correctBoundaryValues(
			[[maybe_unused]] volumeField<tensor3> & corField) const noexcept;

	/*Copy surface values*/
	void correctBoundaryValues(
			[[maybe_unused]] surfaceField<scalar> & corField) const noexcept;

	void correctBoundaryValues(
			[[maybe_unused]] surfaceField<vector> & corField) const noexcept;

	void correctBoundaryValues(
			[[maybe_unused]] surfaceField<tensor> & corField) const noexcept;

	void correctBoundaryValues(
			[[maybe_unused]] surfaceField<tensor3> & corField) const noexcept;

	template<typename T>
	void correctBoundaryConditions(
			[[maybe_unused]] std::array<std::vector<subPatchData<T>>, 6> & conditionsArray) const noexcept
	{
#ifdef MPI_VERSION
		if (!(envNodes.tail < 0))
		{
			auto & conditionsArray_0 = conditionsArray[0];
			for (auto & subp : conditionsArray_0)
				subp.bType = boundaryConditionType::calculatedParallelBoundary;
		}
		if (!(envNodes.point < 0))
		{
			auto & conditionsArray_1 = conditionsArray[1];
			for (auto & subp : conditionsArray_1)
				subp.bType = boundaryConditionType::calculatedParallelBoundary;
		}
		if (!(envNodes.bottom < 0))
		{
			auto & conditionsArray_2 = conditionsArray[2];
			for (auto & subp : conditionsArray_2)
				subp.bType = boundaryConditionType::calculatedParallelBoundary;
		}
		if (!(envNodes.right < 0))
		{
			auto & conditionsArray_3 = conditionsArray[3];
			for (auto & subp : conditionsArray_3)
				subp.bType = boundaryConditionType::calculatedParallelBoundary;
		}
		if (!(envNodes.left < 0))
		{
			auto & conditionsArray_4 = conditionsArray[4];
			for (auto & subp : conditionsArray_4)
				subp.bType = boundaryConditionType::calculatedParallelBoundary;
		}
		if (!(envNodes.top < 0))
		{
			auto & conditionsArray_5 = conditionsArray[5];
			for (auto & subp : conditionsArray_5)
				subp.bType = boundaryConditionType::calculatedParallelBoundary;
		}
#endif
	}

	template<typename typeOfValue>
	void correctBoundaryConditions(
			[[maybe_unused]] volumeField<typeOfValue> & corField) const noexcept
	{
#ifdef MPI_VERSION
		if (!(envNodes.tail < 0))
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.point < 0))
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.bottom < 0))
			for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.right < 0))
			for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.left < 0))
			for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.top < 0))
			for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
#endif
	}

	template<typename typeOfValue>
	void correctBoundaryConditions(
			[[maybe_unused]] surfaceField<typeOfValue> & corField) const noexcept
	{
#ifdef MPI_VERSION
		if (!(envNodes.tail < 0))
			for (std::size_t i = tailSurfaceStart; i < tailSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.point < 0))
			for (std::size_t i = pointSurfaceStart; i < pointSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.bottom < 0))
			for (std::size_t i = bottomSurfaceStart; i < bottomSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.right < 0))
			for (std::size_t i = rightSurfaceStart; i < rightSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.left < 0))
			for (std::size_t i = leftSurfaceStart; i < leftSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
		if (!(envNodes.top < 0))
			for (std::size_t i = topSurfaceStart; i < topSurfaceEnd; ++i)
				corField.boundCond_r()[i].first =
						boundaryConditionType::calculatedParallelBoundary;
#endif
	}

	template<typename typeOfEntiy>
	void correctBoundaryValues(
			[[maybe_unused]] bunchOfFields<typeOfEntiy> & fieldsPack) const noexcept
	{
#ifdef MPI_VERSION
		for (std::size_t k = 1; k < fieldsPack.concentration.v.size(); ++k)
			correctBoundaryValues(fieldsPack.concentration.v[k]);

		correctBoundaryValues(fieldsPack.velocity);
		correctBoundaryValues(fieldsPack.pressure);
		correctBoundaryValues(fieldsPack.kTurb);
		correctBoundaryValues(fieldsPack.epsTurb);
		correctBoundaryValues(fieldsPack.aTurb);
		correctBoundaryValues(fieldsPack.bTurb);
#endif
	}

	void gatherField(const volumeField<scalar> & inField,
			std::valarray<scalar> & retField) noexcept;

	bool isRoot();

	const std::array<std::size_t, 3>& rankDirection() const noexcept;
	const std::array<std::size_t, 3>& rankDirectionMax() const noexcept;

	const bool parallel;
};
}  // namespace schemi

#endif /* MPIHANDLER_HPP_ */
