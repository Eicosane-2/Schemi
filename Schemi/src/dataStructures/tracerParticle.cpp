/*
 * tracerParticle.cpp
 *
 *  Created on: 2025/06/27
 *      Author: Maxim Boldyrev
 */

#include "tracerParticle.hpp"

#include <algorithm>

#include "exception.hpp"
#include "globalConstants.hpp"

schemi::tracerParticle::tracerParticle() :
		position(0), position_1(position), velocity { }
{
}

schemi::tracerParticle::tracerParticle(const vector & inPos,
		const vector & inVelocity) :
		position(inPos), position_1(position), velocity { inVelocity }, step(0)
{
}

schemi::tracerParticle::tracerParticle(const vector & inPos,
		const vector & inPos1, const std::array<vector, 4> & inVelocity,
		const std::size_t inStep) :
		position(inPos), position_1(inPos1), velocity { inVelocity }, step(
				inStep)
{
}

schemi::tracerParticle::~tracerParticle() noexcept
{
}

void schemi::tracerParticle::writeOutput(std::ofstream & output) const
{
	output << std::get<0>(position()) << '\t' << std::get<1>(position()) << '\t'
			<< std::get<2>(position()) << '\n';
	output << std::get<0>(position_1()) << '\t' << std::get<1>(position_1())
			<< '\t' << std::get<2>(position_1()) << '\n';
	for (std::size_t i = 0; i < velocity.size(); ++i)
		output << std::get<0>(velocity[i]()) << '\t'
				<< std::get<1>(velocity[i]()) << '\t'
				<< std::get<2>(velocity[i]()) << '\n';

	output << step << '\n';
}

void schemi::tracerParticle::timeIntegration(const vector & inVelocity,
		const scalar timestep) noexcept
{
	if (step == 0)
	{
		std::get<1>(velocity) = std::get<0>(velocity);
		std::get<0>(velocity) = inVelocity;

		position_1 = position;
		step++;
		position = position_1 + timestep * std::get<0>(velocity);
	}
	else if (step == 1)
	{
		std::get<2>(velocity) = std::get<1>(velocity);
		std::get<1>(velocity) = std::get<0>(velocity);
		std::get<0>(velocity) = inVelocity;

		position_1 = position;
		step++;
		position = position_1
				+ timestep
						* (3. / 2. * std::get<0>(velocity)
								- 1. / 2. * std::get<1>(velocity));
	}
	else if (step == 2)
	{
		std::get<3>(velocity) = std::get<2>(velocity);
		std::get<2>(velocity) = std::get<1>(velocity);
		std::get<1>(velocity) = std::get<0>(velocity);
		std::get<0>(velocity) = inVelocity;

		position_1 = position;
		step++;
		position = position_1
				+ timestep
						* (23. / 12. * std::get<0>(velocity)
								- 16. / 12. * std::get<1>(velocity)
								+ 5. / 12. * std::get<2>(velocity));
	}
	else
	{
		std::get<3>(velocity) = std::get<2>(velocity);
		std::get<2>(velocity) = std::get<1>(velocity);
		std::get<1>(velocity) = std::get<0>(velocity);
		std::get<0>(velocity) = inVelocity;

		position_1 = position;
		position = position_1
				+ timestep
						* (55. / 24. * std::get<0>(velocity)
								- 59. / 24. * std::get<1>(velocity)
								+ 37. / 24. * std::get<2>(velocity)
								- 9. / 24. * std::get<3>(velocity));
	}
}
