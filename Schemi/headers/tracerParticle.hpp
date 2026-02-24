/*
 * tracerParticle.hpp
 *
 *  Created on: 2025/06/27
 *      Author: Maxim Boldyrev
 */

#ifndef TRACERPARTICLE_HPP_
#define TRACERPARTICLE_HPP_

#include <fstream>

#include "mesh.hpp"
#include "vector.hpp"
#include "volumeField.hpp"
#include "surfaceField.hpp"

namespace schemi
{
class tracerParticle
{
	vector position { 0 }, position_1 { 0 };
	std::array<vector, 4> velocity { };
	std::size_t step { 0 };

protected:
	virtual void writeOutput(std::ofstream & output) const;
	void timeIntegration(const vector & inVelocity,
			const scalar timestep) noexcept;
public:
	tracerParticle() noexcept = default;
	tracerParticle(const vector & inPos, const vector & inVelocity) noexcept;
	tracerParticle(const vector & inPos, const vector & inPos1,
			const std::array<vector, 4> & inVelocity,
			const std::size_t inStep) noexcept;
	virtual ~tracerParticle() noexcept = default;

	const vector& getPosition() const noexcept
	{
		return position;
	}
	const std::array<vector, 4>& getVelocities() const noexcept
	{
		return velocity;
	}
};
}  // namespace schemi

#endif /* TRACERPARTICLE_HPP_ */
