/*
 * GoncharovTracerModel.cpp
 *
 *  Created on: 2025/06/28
 *      Author: Maxim Boldyrev
 */

#include "GoncharovTracerModel.hpp"

#include "globalConstants.hpp"
#include "vector.hpp"
#include "intExpPow.hpp"

void schemi::GoncharovTracerModel::initialiseTransitionCheck(
		const std::string & initCheckMethod)
{
	if (initCheckMethod == "viaRe")
		transitioniCriterion = [this](const scalar Re_in, const scalar) -> bool 
		{	return Re_in>=this->ReCriterion;};
	else if (initCheckMethod == "viaEta")
		transitioniCriterion = [this](const scalar,
				const scalar dEtaL_in) -> bool 
				{	return dEtaL_in>=this->relativeCriterion;};
	else if (initCheckMethod == "viaReOrEta")
		transitioniCriterion =
				[this](const scalar Re_in,
						const scalar dEtaL_in) -> bool 
						{	return (Re_in>=this->ReCriterion) || (dEtaL_in>=this->relativeCriterion);};
	else if (initCheckMethod == "viaReAndEta")
		transitioniCriterion =
				[this](const scalar Re_in,
						const scalar dEtaL_in) -> bool 
						{	return (Re_in>=this->ReCriterion) && (dEtaL_in>=this->relativeCriterion);};
	else
		[[unlikely]]
		throw exception("unknown type of turbulence initialisation criterion.",
				errors::initialisationError);
}

schemi::GoncharovTracerModel::GoncharovTracerModel() noexcept :
		tracerParticle()
{
}

schemi::GoncharovTracerModel::GoncharovTracerModel(
		const std::string & initCheckMethod, const vector & inPos,
		const vector & inVelocity, const std::size_t sub1,
		const std::size_t sub2, const int pertType, const scalar eta0In,
		const scalar lambdaIn, const scalar radiusOfIfluenceIn) :
		tracerParticle(inPos, inVelocity), s12 { sub1, sub2 }, c(pertType), eta_0(
				eta0In), lambda(lambdaIn), k(2 * Pi_number / lambda), radiusOfIfluence(
				radiusOfIfluenceIn), eta(eta_0), eta_1(eta_0), eta_2(eta_0)
{
	initialiseTransitionCheck(initCheckMethod);
}

schemi::GoncharovTracerModel::GoncharovTracerModel(
		const std::string & initCheckMethod, const vector & inPos,
		const vector & inPos1, const std::array<vector, 4> & inVelocity,
		const std::size_t inStep, const std::size_t sub1,
		const std::size_t sub2, const int pertType, const scalar eta0In,
		const scalar lambdaIn, const scalar radiusOfIfluenceIn,
		const scalar etaCur, const scalar eta1Cur, const scalar eta2Cur,
		const scalar rho1Cur, const scalar rho2Cur, const scalar k0Cur,
		const scalar eps0Cur, const scalar b0Cur, const vector & a0Cur,
		const interfaceStatus curStatus) :
		tracerParticle(inPos, inPos1, inVelocity, inStep), s12 { sub1, sub2 }, c(
				pertType), eta_0(eta0In), lambda(lambdaIn), k(
				2 * Pi_number / lambda), radiusOfIfluence(radiusOfIfluenceIn), eta(
				etaCur), eta_1(eta1Cur), eta_2(eta2Cur), rho1(rho1Cur), rho2(
				rho2Cur), k0(k0Cur), eps0(eps0Cur), b0(b0Cur), a0(a0Cur), status(
				curStatus)
{
	initialiseTransitionCheck(initCheckMethod);
}

schemi::GoncharovTracerModel::~GoncharovTracerModel() noexcept
{
}

void schemi::GoncharovTracerModel::writeOutput(std::ofstream & output) const
{
	tracerParticle::writeOutput(output);

	output << eta << '\t' << eta_1 << '\t' << eta_2 << '\n';
	output << rho1 << '\t' << rho2 << '\n';
	output << k0 << '\t' << eps0 << '\t' << b0 << '\n';
	output << std::get<0>(a0()) << '\t' << std::get<1>(a0()) << '\t'
			<< std::get<2>(a0()) << '\n';
	output << static_cast<int>(status) << '\n';
}

void schemi::GoncharovTracerModel::timeIntegration(const scalar density1,
		const scalar density2, const vector & gradRho, const vector & u,
		const scalar timestep, const vector&, const scalar, const scalar,
		const tensor&)
{
	tracerParticle::timeIntegration(u, timestep);

	if (status == interfaceStatus::notDeveloped)
	{
		rho1 = density1;
		rho2 = density2;

		eta_2 = eta_1;
		eta_1 = eta;

		const auto normale = gradRho / -gradRho.mag(); // From heavy to light.

		const scalar At = std::abs(rho2 - rho1) / (rho1 + rho2);

		const auto gVec = (std::get<0>(getVelocities())
				- std::get<1>(getVelocities())) / timestep;

		const auto eta2 = -c * k / (4 * (1 + c))
				* (1
						+ ((1 + c) * eta_0 * k - 1)
								* std::exp(-k * (1 + c) * (eta_1 - eta_0)));

		const auto F1 = 2 * At * pow<scalar, 2>(eta2)
				+ pow<scalar, 2>(c) * At * k * eta2 / (2 * (1 + c))
				- pow<scalar, 2>(c * k) / (8 * (1 + c));

		const auto F2 = 2 * At * pow<scalar, 2>(eta2)
				+ (At + c * At - 2 * c - 1) * k * eta2 / (1 + c)
				+ c * pow<scalar, 2>(k) * (3 * c * At / 2 + At - c - 1)
						/ (4 * pow<scalar, 2>(1 + c));

		const auto D = eta2 - c * k / (4 * (1 + c));

		const auto g = std::abs(gVec & normale);

		const auto A = F1 / D;
		const auto B = F2 * pow<scalar, 2>(c * k / D) / 8 / A;
		const auto C = 2 * g * At * eta2 / A;

		eta = 2 * eta_1 - eta_2 - pow<scalar, 2>(eta_1 - eta_2) * B
				- C * pow<scalar, 2>(timestep);

		if (eta < eta_0)
			throw exception("Perturbation amplitude less than initial.",
					errors::systemError);
	}
}

schemi::interfaceStatus schemi::GoncharovTracerModel::checkTransition(
		const scalar nu, const scalar timestep, const vector & cellRadius,
		const vector & surfaceRadius) noexcept
{
	if (status != interfaceStatus::notDeveloped)
		return status;

	const auto detadt = std::abs(eta - eta_1) / timestep;

	Re = eta / nu * detadt;
	relDeltaEta = std::abs(eta - eta_0) / lambda;

	if (transitioniCriterion(Re, relDeltaEta))
	{
		const auto charactDist = (cellRadius - surfaceRadius).mag();

		if ((cellCoefficient * charactDist) > eta)
			status = interfaceStatus::developedNotResolvable;
		else
			status = interfaceStatus::developedResolvable;
	}
	else
		status = interfaceStatus::notDeveloped;

	if (status == interfaceStatus::notDeveloped)
		return status;
	else
	{
		const auto denWeighted = rho1 * rho2 / pow<scalar, 2>(rho1 + rho2);
		const auto deltaDen = pow<scalar, 2>(rho2 - rho1) / (rho1 * rho2);

		k0 = Ck * 1.5 * denWeighted
				* pow<scalar, 2>(std::get<0>(getVelocities()).mag());

		const auto s0 = 2 * Ceps * eta;
		eps0 = pow<scalar, 3>(std::sqrt(k0)) / s0;

		b0 = Cb / 4 * deltaDen;

		return status;
	}
}
