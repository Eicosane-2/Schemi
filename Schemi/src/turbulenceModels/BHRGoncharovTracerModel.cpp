/*
 * BHRGoncharovTracerModel.cpp
 *
 *  Created on: 2025/12/06
 *      Author: Maxim Boldyrev
 */

#include "BHRGoncharovTracerModel.hpp"

#include "intExpPow.hpp"

schemi::BHRGoncharovTracerModel::BHRGoncharovTracerModel() noexcept :
		GoncharovTracerModel()
{
}

schemi::BHRGoncharovTracerModel::BHRGoncharovTracerModel(
		const std::string & initCheckMethod, const vector & inPos,
		const vector & inVelocity, const std::size_t sub1,
		const std::size_t sub2, const int pertType, const scalar eta0In,
		const scalar lambdaIn, const scalar radiusOfIfluenceIn,
		const scalar CmuIn, const scalar C0In, const scalar C2In,
		const scalar C3In, const scalar Ca1In, const scalar Cb1In) :
		GoncharovTracerModel(initCheckMethod, inPos, inVelocity, sub1, sub2,
				pertType, eta0In, lambdaIn, radiusOfIfluenceIn), Cmu(CmuIn), C0(
				C0In), C2(C2In), C3(C3In), Ca1(Ca1In), Cb1(Cb1In)
{
}

schemi::BHRGoncharovTracerModel::BHRGoncharovTracerModel(
		const std::string & initCheckMethod, const vector & inPos,
		const vector & inPos1, const std::array<vector, 4> & inVelocity,
		const std::size_t inStep, const std::size_t sub1,
		const std::size_t sub2, const int pertType, const scalar eta0In,
		const scalar lambdaIn, const scalar radiusOfIfluenceIn,
		const scalar etaCur, const scalar eta1Cur, const scalar eta2Cur,
		const scalar rho1Cur, const scalar rho2Cur, const scalar k0Cur,
		const scalar eps0Cur, const scalar b0Cur, const vector & a0Cur,
		const interfaceStatus curStatus, const scalar CmuIn, const scalar C0In,
		const scalar C2In, const scalar C3In, const scalar Ca1In,
		const scalar Cb1In) :
		GoncharovTracerModel(initCheckMethod, inPos, inPos1, inVelocity, inStep,
				sub1, sub2, pertType, eta0In, lambdaIn, radiusOfIfluenceIn,
				etaCur, eta1Cur, eta2Cur, rho1Cur, rho2Cur, k0Cur, eps0Cur,
				b0Cur, a0Cur, curStatus), Cmu(CmuIn), C0(C0In), C2(C2In), C3(
				C3In), Ca1(Ca1In), Cb1(Cb1In)
{
}

schemi::BHRGoncharovTracerModel::~BHRGoncharovTracerModel() noexcept
{
}

void schemi::BHRGoncharovTracerModel::timeIntegration(const scalar density1,
		const scalar density2, const vector & gradRho, const vector & u,
		const vector & g, const scalar timestep, const vector & gradP,
		const scalar rho, const scalar divU, const tensor & gradU)
{
	const auto currentStatus { getStatus() };

	if (currentStatus == interfaceStatus::notDeveloped)
	{
		GoncharovTracerModel::timeIntegration(density1, density2, gradRho, u, g,
				timestep);
	}
	else if (currentStatus == interfaceStatus::developedNotResolvable)
	{
		tracerParticle::timeIntegration(u, timestep);

		const auto kOld = getk0();
		const auto epsOld = geteps0();
		const auto aOld = geta0();
		const auto bOld = getb0();
		const auto etaOld = getEta();

		const vector gradP_rho(gradP / rho);
		const vector gradRho_rho(gradRho / rho);
		const scalar aGradPRho = aOld & gradP_rho;
		const scalar kDivU = twothirds * kOld * divU;
		const scalar nut = Cmu * pow<scalar, 2>(kOld) / epsOld;
		const scalar diffusion = nut * pow<scalar, 2>(Ceta / etaOld);
		const scalar ek = epsOld / kOld;

		const scalar kSource = aGradPRho - epsOld - kDivU - kOld * diffusion;

		const scalar epsSource = ek
				* (C0 * aGradPRho - C2 * epsOld - C3 * kDivU)
				- epsOld * diffusion;

		const vector aSource = bOld * gradP_rho - kOld * gradRho_rho
				- (aOld & gradU) - aOld * diffusion - Ca1 * ek * aOld;

		const scalar bSource = -2 * (bOld + 1) * (aOld & gradRho_rho)
				- bOld * diffusion - Cb1 * ek * bOld;

		const scalar etaSource = etaOld * diffusion + etaOld * divU;

		const auto kNew = kOld + kSource * timestep;
		const auto epsNew = epsOld + epsSource * timestep;
		const auto aNew = aOld + aSource * timestep;
		const auto bNew = bOld + bSource * timestep;
		const auto etaNew = etaOld + etaSource * timestep;

		setEta(etaNew);
		setk0(kNew);
		seteps0(epsNew);
		seta0(aNew);
		setb0(bNew);
	}
	else
		tracerParticle::timeIntegration(u, timestep);
}

schemi::interfaceStatus schemi::BHRGoncharovTracerModel::checkTransition(
		const scalar nu, const scalar timestep, const vector & cellRadius,
		const vector & surfaceRadius) noexcept
{
	auto currentstatus = getStatus();

	switch (currentstatus)
	{
	case interfaceStatus::developedResolvable:
		break;
	case interfaceStatus::notDeveloped:
		currentstatus = GoncharovTracerModel::checkTransition(nu, timestep,
				cellRadius, surfaceRadius);
		break;
	case interfaceStatus::developedNotResolvable:
	default:
	{
		const auto charactDist = (cellRadius - surfaceRadius).mag();

		if ((cellCoefficient * charactDist) > getEta())
			currentstatus = interfaceStatus::developedNotResolvable;
		else
			currentstatus = interfaceStatus::developedResolvable;

		setStatus(currentstatus);
	}
		break;
	}

	return currentstatus;
}
