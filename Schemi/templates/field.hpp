/*
 * field.hpp
 *
 *  Created on: 2019/11/13
 *      Author: Maxim Boldyrev
 *
 *      Structure storing physical fields.
 */
#ifndef FIELD_HPP_
#define FIELD_HPP_

#include <valarray>

#include "exception.hpp"
#include "globalConstants.hpp"
#include "mesh.hpp"
#include "subPatchData.hpp"

namespace schemi
{
template<typename typeOfValue, typename typeOfEntity>
struct field
{
	field(const field<typeOfValue, typeOfEntity>&) = default;

	auto& operator=(const field<typeOfValue, typeOfEntity> & f)
	{
		if ((f.boundCond().size() != boundaryConditionInfo.size())
				&& (boundaryConditionInfo.size() != 0))
			throw exception("Fields have different boundary condition size.",
					errorsEnum::fieldInitializationError);

		valueField = f.ref();
		boundaryConditionInfo = f.boundCond();

		fieldSize = f.size();

		return *this;
	}

	explicit field(const mesh & meshIn,
			const typeOfValue & value = typeOfValue(0)) :
			valueField(0), boundaryConditionInfo(0), meshReference { meshIn }, fieldSize {
					0 }
	{
		is_mesh_initialised();

		if constexpr (std::is_same_v<typeOfEntity, cubicCell>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.cellsSize());

			boundaryConditionInfo.resize(meshReference.surfacesSize());
			for (std::size_t i = 0; i < boundaryConditionInfo.size(); ++i)
			{
				boundaryConditionInfo[i].first = meshReference.bndType()[i];
				boundaryConditionInfo[i].second = value;
			}

			fieldSize = meshReference.cellsSize();
		}
		else if constexpr (std::is_same_v<typeOfEntity, quadraticSurface>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.surfacesSize());

			fieldSize = meshReference.surfacesSize();
		}
		else
			throw exception("Unknown type of field.",
					errorsEnum::fieldInitializationError);
	}

	field(const mesh & meshIn, const typeOfValue & value,
			const subPatchData<typeOfValue> & bData1,
			const subPatchData<typeOfValue> & bData2,
			const subPatchData<typeOfValue> & bData3,
			const subPatchData<typeOfValue> & bData4,
			const subPatchData<typeOfValue> & bData5,
			const subPatchData<typeOfValue> & bData6) :
			valueField(0), boundaryConditionInfo(0), meshReference { meshIn }, fieldSize {
					0 }
	{
		is_mesh_initialised();

		if constexpr (std::is_same_v<typeOfEntity, cubicCell>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.cellsSize());

			const std::size_t surfSize = meshReference.surfacesSize();
			boundaryConditionInfo.resize(surfSize,
					std::pair<boundaryConditionType, typeOfValue>(
							boundaryConditionType::innerSurface,
							typeOfValue(0)));

			std::size_t prev { 0 };

			for (std::size_t j = prev; j < (prev + meshReference.tailNumber());
					++j)
			{
				boundaryConditionInfo[j].first = bData1.bType;
				boundaryConditionInfo[j].second = bData1.fixVal;
			}
			prev += meshReference.tailNumber() + meshReference.innerNumber();

			for (std::size_t j = prev; j < (prev + meshReference.pointNumber());
					++j)
			{
				boundaryConditionInfo[j].first = bData2.bType;
				boundaryConditionInfo[j].second = bData2.fixVal;
			}
			prev += meshReference.pointNumber();

			for (std::size_t j = prev;
					j < (prev + meshReference.bottomNumber()); ++j)
			{
				boundaryConditionInfo[j].first = bData3.bType;
				boundaryConditionInfo[j].second = bData3.fixVal;
			}
			prev += meshReference.bottomNumber();

			for (std::size_t j = prev; j < (prev + meshReference.rightNumber());
					++j)
			{
				boundaryConditionInfo[j].first = bData4.bType;
				boundaryConditionInfo[j].second = bData4.fixVal;
			}
			prev += meshReference.rightNumber();

			for (std::size_t j = prev; j < (prev + meshReference.leftNumber());
					++j)
			{
				boundaryConditionInfo[j].first = bData5.bType;
				boundaryConditionInfo[j].second = bData5.fixVal;
			}
			prev += meshReference.leftNumber();

			for (std::size_t j = prev; j < (prev + meshReference.topNumber());
					++j)
			{
				boundaryConditionInfo[j].first = bData6.bType;
				boundaryConditionInfo[j].second = bData6.fixVal;
			}

			fieldSize = meshReference.cellsSize();
		}
		else
			staticNotVolumeField();
	}

	field(const mesh & meshIn, const typeOfValue & value,
			const std::vector<subPatchData<typeOfValue>> & bData1,
			const std::vector<subPatchData<typeOfValue>> & bData2,
			const std::vector<subPatchData<typeOfValue>> & bData3,
			const std::vector<subPatchData<typeOfValue>> & bData4,
			const std::vector<subPatchData<typeOfValue>> & bData5,
			const std::vector<subPatchData<typeOfValue>> & bData6) :
			valueField(0), boundaryConditionInfo(0), meshReference { meshIn }, fieldSize {
					0 }
	{
		is_mesh_initialised();

		if constexpr (std::is_same_v<typeOfEntity, cubicCell>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.cellsSize());

			const std::size_t surfSize = meshReference.surfacesSize();
			boundaryConditionInfo.resize(surfSize,
					std::pair<boundaryConditionType, typeOfValue>(
							boundaryConditionType::innerSurface,
							typeOfValue(0)));

			std::size_t prev { 0 };

			for (std::size_t j = prev; j < (prev + meshReference.tailNumber());
					++j)
			{
				const vector surfaceCoord = meshReference.surfaces()[j].rC();

				bool zoneFounded { false };

				for (std::size_t zone = 0; zone < bData1.size(); ++zone)
				{
					if (true
							&& (surfaceCoord.v()[1]
									>= bData1[zone].patchBeg.v()[1])
							&& (surfaceCoord.v()[1]
									<= bData1[zone].patchEnd.v()[1])
							&& (surfaceCoord.v()[2]
									>= bData1[zone].patchBeg.v()[2])
							&& (surfaceCoord.v()[2]
									<= bData1[zone].patchEnd.v()[2]))
					{
						zoneFounded = true;
						boundaryConditionInfo[j].first = bData1[zone].bType;
						boundaryConditionInfo[j].second = bData1[zone].fixVal;
						break;
					}
				}

				if (!zoneFounded)
					throw exception(
							"Surface " + std::to_string(j)
									+ " of a dropped out of all zones.",
							errorsEnum::initializationError);
			}
			prev += meshReference.tailNumber() + meshReference.innerNumber();

			for (std::size_t j = prev; j < (prev + meshReference.pointNumber());
					++j)
			{
				const vector surfaceCoord = meshReference.surfaces()[j].rC();

				bool zoneFounded { false };

				for (std::size_t zone = 0; zone < bData2.size(); ++zone)
				{
					if (true
							&& (surfaceCoord.v()[1]
									>= bData2[zone].patchBeg.v()[1])
							&& (surfaceCoord.v()[1]
									<= bData2[zone].patchEnd.v()[1])
							&& (surfaceCoord.v()[2]
									>= bData2[zone].patchBeg.v()[2])
							&& (surfaceCoord.v()[2]
									<= bData2[zone].patchEnd.v()[2]))
					{
						zoneFounded = true;
						boundaryConditionInfo[j].first = bData2[zone].bType;
						boundaryConditionInfo[j].second = bData2[zone].fixVal;
						break;
					}
				}

				if (!zoneFounded)
					throw exception(
							"Surface " + std::to_string(j)
									+ " of a dropped out of all zones.",
							errorsEnum::initializationError);
			}
			prev += meshReference.pointNumber();

			for (std::size_t j = prev;
					j < (prev + meshReference.bottomNumber()); ++j)
			{
				const vector surfaceCoord = meshReference.surfaces()[j].rC();

				bool zoneFounded { false };

				for (std::size_t zone = 0; zone < bData3.size(); ++zone)
				{
					if ((surfaceCoord.v()[0] >= bData3[zone].patchBeg.v()[0])
							&& (surfaceCoord.v()[0]
									<= bData3[zone].patchEnd.v()[0])
							&& (surfaceCoord.v()[1]
									>= bData3[zone].patchBeg.v()[1])
							&& (surfaceCoord.v()[1]
									<= bData3[zone].patchEnd.v()[1]) && true)
					{
						zoneFounded = true;
						boundaryConditionInfo[j].first = bData3[zone].bType;
						boundaryConditionInfo[j].second = bData3[zone].fixVal;
						break;
					}
				}

				if (!zoneFounded)
					throw exception(
							"Surface " + std::to_string(j)
									+ " of a dropped out of all zones.",
							errorsEnum::initializationError);
			}
			prev += meshReference.bottomNumber();

			for (std::size_t j = prev; j < (prev + meshReference.rightNumber());
					++j)
			{
				const vector surfaceCoord = meshReference.surfaces()[j].rC();

				bool zoneFounded { false };

				for (std::size_t zone = 0; zone < bData4.size(); ++zone)
				{
					if ((surfaceCoord.v()[0] >= bData4[zone].patchBeg.v()[0])
							&& (surfaceCoord.v()[0]
									<= bData4[zone].patchEnd.v()[0]) && true
							&& (surfaceCoord.v()[2]
									>= bData4[zone].patchBeg.v()[2])
							&& (surfaceCoord.v()[2]
									<= bData4[zone].patchEnd.v()[2]))
					{
						zoneFounded = true;
						boundaryConditionInfo[j].first = bData4[zone].bType;
						boundaryConditionInfo[j].second = bData4[zone].fixVal;
						break;
					}
				}

				if (!zoneFounded)
					throw exception(
							"Surface " + std::to_string(j)
									+ " of a dropped out of all zones.",
							errorsEnum::initializationError);
			}
			prev += meshReference.rightNumber();

			for (std::size_t j = prev; j < (prev + meshReference.leftNumber());
					++j)
			{
				const vector surfaceCoord = meshReference.surfaces()[j].rC();

				bool zoneFounded { false };

				for (std::size_t zone = 0; zone < bData5.size(); ++zone)
				{
					if ((surfaceCoord.v()[0] >= bData5[zone].patchBeg.v()[0])
							&& (surfaceCoord.v()[0]
									<= bData5[zone].patchEnd.v()[0]) && true
							&& (surfaceCoord.v()[2]
									>= bData5[zone].patchBeg.v()[2])
							&& (surfaceCoord.v()[2]
									<= bData5[zone].patchEnd.v()[2]))
					{
						zoneFounded = true;
						boundaryConditionInfo[j].first = bData5[zone].bType;
						boundaryConditionInfo[j].second = bData5[zone].fixVal;
						break;
					}
				}

				if (!zoneFounded)
					throw exception(
							"Surface " + std::to_string(j)
									+ " of a dropped out of all zones.",
							errorsEnum::initializationError);
			}
			prev += meshReference.leftNumber();

			for (std::size_t j = prev; j < (prev + meshReference.topNumber());
					++j)
			{
				const vector surfaceCoord = meshReference.surfaces()[j].rC();

				bool zoneFounded { false };

				for (std::size_t zone = 0; zone < bData6.size(); ++zone)
				{
					if ((surfaceCoord.v()[0] >= bData6[zone].patchBeg.v()[0])
							&& (surfaceCoord.v()[0]
									<= bData6[zone].patchEnd.v()[0])
							&& (surfaceCoord.v()[1]
									>= bData6[zone].patchBeg.v()[1])
							&& (surfaceCoord.v()[1]
									<= bData6[zone].patchEnd.v()[1]) && true)
					{
						zoneFounded = true;
						boundaryConditionInfo[j].first = bData6[zone].bType;
						boundaryConditionInfo[j].second = bData6[zone].fixVal;
						break;
					}
				}

				if (!zoneFounded)
					throw exception(
							"Surface " + std::to_string(j)
									+ " of a dropped out of all zones.",
							errorsEnum::initializationError);
			}

			fieldSize = meshReference.cellsSize();
		}
		else
			staticNotVolumeField();
	}

	field(const mesh & meshIn, const typeOfValue & value,
			const std::vector<std::pair<boundaryConditionType, typeOfValue>> & boundCondtIn) :
			valueField(0), boundaryConditionInfo(0), meshReference { meshIn }, fieldSize {
					0 }
	{
		is_mesh_initialised();

		if constexpr (std::is_same_v<typeOfEntity, cubicCell>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.cellsSize());

			fieldSize = meshReference.cellsSize();
		}
		else if constexpr (std::is_same_v<typeOfEntity, quadraticSurface>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.surfacesSize());

			fieldSize = meshReference.surfacesSize();
		}
		else
			throw exception("Unknown type of field.",
					errorsEnum::fieldInitializationError);

		boundaryConditionInfo = boundCondtIn;
	}

	const std::valarray<typeOfValue>& ref() const noexcept
	{
		return valueField;
	}

	std::valarray<typeOfValue>& ref_r() noexcept
	{
		return valueField;
	}

	const std::vector<std::pair<boundaryConditionType, typeOfValue>>& boundCond() const noexcept
	{
		return boundaryConditionInfo;
	}

	std::vector<std::pair<boundaryConditionType, typeOfValue>>& boundCond_r() noexcept
	{
		return boundaryConditionInfo;
	}

	const mesh& meshRef() const noexcept
	{
		return meshReference;
	}

	std::size_t size() const noexcept
	{
		return fieldSize;
	}
private:
	std::valarray<typeOfValue> valueField;
	std::vector<std::pair<boundaryConditionType, typeOfValue>> boundaryConditionInfo;
	const mesh & meshReference;
	std::size_t fieldSize;

	void is_mesh_initialised() const
	{
		if (!meshReference.is_initialized())
			throw exception("Mesh has not yet been initialized.",
					errorsEnum::fieldInitializationError);
	}

	template<bool flag = false>
	constexpr void staticNotVolumeField()
	{
		static_assert(flag, "Field is not volumeField.");
	}
};
}  // namespace schemi

#endif /* FIELD_HPP_ */
