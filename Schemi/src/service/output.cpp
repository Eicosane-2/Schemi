/*
 * output.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "output.hpp"

#include <fstream>
#include <iostream>
#include <filesystem>

void schemi::output::dataOutput(const structForOutput & outputData,
		const std::size_t noutput, const scalar Time)
{
	const auto noutputStr = std::to_string(noutput);
	const std::size_t length_of_noutput = noutputStr.size();

	if (lengthOfNumber < length_of_noutput)
		throw exception(
				"Output number length larger than specified number of digits.",
				errors::tooBigOutputNumberError);

	std::string bufOutputN(lengthOfNumber - length_of_noutput, '0');

	bufOutputN.append(noutputStr);

	std::string outputFileName { "./result/output_" };

	std::string timeFileName { "./result/Time.tsv" };

	outputFileName.append(bufOutputN);
	outputFileName.append(".tsv");

	std::ofstream outputFile { outputFileName };
	outputFile.precision(ioPrecision);

	std::ofstream timeFile { timeFileName, std::ios::app };
	timeFile.precision(ioPrecision);

	if (outputFile.is_open())
		std::cout << outputFileName
				<< " output file is opened. Number of output: " << noutput
				<< '.' << std::endl;
	else
		[[unlikely]]
		throw std::ofstream::failure(
				std::string("Couldn't create outputFile ")
						+ std::string(outputFileName) + std::string("."));

	if (timeFile.is_open())
		std::cout << "./result/Time.tsv is opened." << std::endl;
	else
		[[unlikely]]
		throw std::ofstream::failure("Couldn't open ./result/Time.tsv.");

	outputFile << "coordinate_x" << '\t';
	outputFile << "coordinate_y" << '\t';
	outputFile << "coordinate_z" << '\t';
	for (std::size_t k = 1; k < outputData.concentration.size(); ++k)
	{
		std::string concentrationName { "concentration_" }, bufComponentO(
				std::to_string(k));

		concentrationName.append(bufComponentO);
		outputFile << concentrationName << '\t';
	}

	outputFile << "concentration" << '\t';
	for (std::size_t k = 1; k < outputData.concentration.size(); ++k)
	{
		std::string molarFactionName { "x_" }, bufComponentO(std::to_string(k));

		molarFactionName.append(bufComponentO);
		outputFile << molarFactionName << '\t';
	}

	for (std::size_t k = 1; k < outputData.density.size(); ++k)
	{
		std::string densityName { "density_" }, bufComponentO(
				std::to_string(k));

		densityName.append(bufComponentO);
		outputFile << densityName << '\t';
	}

	outputFile << "density" << '\t';
	for (std::size_t k = 1; k < outputData.density.size(); ++k)
	{
		std::string massFractionName { "w_" }, bufComponentO(std::to_string(k));

		massFractionName.append(bufComponentO);
		outputFile << massFractionName << '\t';
	}

	outputFile << "velocity_x" << '\t';
	outputFile << "velocity_y" << '\t';
	outputFile << "velocity_z" << '\t';
	outputFile << "pressure" << '\t';
	outputFile << "temperature" << '\t';
	outputFile << "k" << '\t';
	outputFile << "epsilon" << '\t';
	outputFile << "a_x" << '\t';
	outputFile << "a_y" << '\t';
	outputFile << "a_z" << '\t';
	outputFile << "b" << '\t';
	outputFile << "turbulent_viscosity" << '\t';
	outputFile << "momentum_x" << '\t';
	outputFile << "momentum_y" << '\t';
	outputFile << "momentum_z" << '\t';
	outputFile << "total_energy" << '\t';
	outputFile << "internal_energy" << '\t';
	outputFile << "Helmholtz_energy" << '\t';
	outputFile << "entropy" << '\t';
	outputFile << "sonic_speed" << '\n';

	for (std::size_t i = 0; i < outputData.x_coord.size(); ++i)
	{
		outputFile << outputData.x_coord[i] << '\t';
		outputFile << outputData.y_coord[i] << '\t';
		outputFile << outputData.z_coord[i] << '\t';

		for (std::size_t k = 1; k < outputData.concentration.size(); ++k)
			outputFile << outputData.concentration[k][i] << '\t';

		outputFile << outputData.concentration[0][i] << '\t';

		for (std::size_t k = 1; k < outputData.concentration.size(); ++k)
			outputFile
					<< outputData.concentration[k][i]
							/ outputData.concentration[0][i] << '\t';

		for (std::size_t k = 1; k < outputData.density.size(); ++k)
			outputFile << outputData.density[k][i] << '\t';

		outputFile << outputData.density[0][i] << '\t';

		for (std::size_t k = 1; k < outputData.density.size(); ++k)
			outputFile << outputData.density[k][i] / outputData.density[0][i]
					<< '\t';

		outputFile << outputData.velocity_x[i] << '\t';
		outputFile << outputData.velocity_y[i] << '\t';
		outputFile << outputData.velocity_z[i] << '\t';

		outputFile << outputData.pressure[i] << '\t';

		outputFile << outputData.temperature[i] << '\t';

		outputFile << outputData.kTurb[i] << '\t';

		outputFile << outputData.epsTurb[i] << '\t';

		outputFile << outputData.aTurb_x[i] << '\t';
		outputFile << outputData.aTurb_y[i] << '\t';
		outputFile << outputData.aTurb_z[i] << '\t';

		outputFile << outputData.bTurb[i] << '\t';

		outputFile << outputData.tNu[i] << '\t';

		outputFile << outputData.momentum_x[i] << '\t';
		outputFile << outputData.momentum_y[i] << '\t';
		outputFile << outputData.momentum_z[i] << '\t';

		outputFile << outputData.totalEnergy[i] << '\t';

		outputFile << outputData.internalEnergy[i] << '\t';

		outputFile << outputData.HelmholtzEnergy[i] << '\t';

		outputFile << outputData.entropy[i] << '\t';

		outputFile << outputData.sonicSpeed[i] << '\n';
	}

	timeFile << noutput << '\t' << Time << '\n';

	outputFile.close();

	timeFile.close();

	std::string fieldDataDirectoryName { "./fieldsOutput/output_" };
	fieldDataDirectoryName.append(bufOutputN);
	fieldDataDirectoryName.append("/");
	const auto fieldDataFileName_a = fieldDataDirectoryName
			+ std::string("a.dat");
	const auto fieldDataFileName_b = fieldDataDirectoryName
			+ std::string("b.dat");
	const auto fieldDataFileName_epsilon = fieldDataDirectoryName
			+ std::string("epsilon.dat");
	const auto fieldDataFileName_k = fieldDataDirectoryName
			+ std::string("k.dat");
	const auto fieldDataFileName_pressure = fieldDataDirectoryName
			+ std::string("pressure.dat");
	const auto fieldDataFileName_velocity = fieldDataDirectoryName
			+ std::string("velocity.dat");
	std::vector<std::string> fieldDataFileName_concentration(
			outputData.concentrationNonSorted.size());
	for (std::size_t k = 0; k < fieldDataFileName_concentration.size(); ++k)
		fieldDataFileName_concentration[k] = fieldDataDirectoryName
				+ std::string("sub_") + std::string(std::to_string(k + 1))
				+ std::string(".dat");

	std::filesystem::create_directories(fieldDataDirectoryName);

	std::ofstream output_aVector { fieldDataFileName_a };
	output_aVector.precision(ioPrecision);
	if (!output_aVector.is_open())
		throw std::ofstream::failure(
				std::string("Couldn't create output file for field data ")
						+ std::string(fieldDataFileName_a) + std::string("."));

	for (std::size_t i = 0; i < outputData.aTurb_xNonSorted.size(); ++i)
		output_aVector << outputData.aTurb_xNonSorted[i] << '\t'
				<< outputData.aTurb_yNonSorted[i] << '\t'
				<< outputData.aTurb_zNonSorted[i] << '\n';

	output_aVector.close();

	std::ofstream output_bScalar { fieldDataFileName_b };
	output_bScalar.precision(ioPrecision);
	if (!output_bScalar.is_open())
		throw std::ofstream::failure(
				std::string("Couldn't create output file for field data ")
						+ std::string(fieldDataFileName_b) + std::string("."));

	for (std::size_t i = 0; i < outputData.bTurbNonSorted.size(); ++i)
		output_bScalar << outputData.bTurbNonSorted[i] << '\n';

	output_bScalar.close();

	std::ofstream output_epsilonScalar { fieldDataFileName_epsilon };
	output_epsilonScalar.precision(ioPrecision);
	if (!output_epsilonScalar.is_open())
		throw std::ofstream::failure(
				std::string("Couldn't create output file for field data ")
						+ std::string(fieldDataFileName_epsilon)
						+ std::string("."));

	for (std::size_t i = 0; i < outputData.epsTurbNonSorted.size(); ++i)
		output_epsilonScalar << outputData.epsTurbNonSorted[i] << '\n';

	output_epsilonScalar.close();

	std::ofstream output_kScalar { fieldDataFileName_k };
	output_kScalar.precision(ioPrecision);
	if (!output_kScalar.is_open())
		throw std::ofstream::failure(
				std::string("Couldn't create output file for field data ")
						+ std::string(fieldDataFileName_k) + std::string("."));

	for (std::size_t i = 0; i < outputData.kTurbNonSorted.size(); ++i)
		output_kScalar << outputData.kTurbNonSorted[i] << '\n';

	output_kScalar.close();

	std::ofstream output_pressureScalar { fieldDataFileName_pressure };
	output_pressureScalar.precision(ioPrecision);
	if (!output_pressureScalar.is_open())
		throw std::ofstream::failure(
				std::string("Couldn't create output file for field data ")
						+ std::string(fieldDataFileName_pressure)
						+ std::string("."));

	for (std::size_t i = 0; i < outputData.pressureNonSorted.size(); ++i)
		output_pressureScalar << outputData.pressureNonSorted[i] << '\n';

	output_pressureScalar.close();

	std::ofstream output_velocityVector { fieldDataFileName_velocity };
	output_velocityVector.precision(ioPrecision);
	if (!output_velocityVector.is_open())
		throw std::ofstream::failure(
				std::string("Couldn't create output file for field data ")
						+ std::string(fieldDataFileName_velocity)
						+ std::string("."));

	for (std::size_t i = 0; i < outputData.velocity_xNonSorted.size(); ++i)
		output_velocityVector << outputData.velocity_xNonSorted[i] << '\t'
				<< outputData.velocity_yNonSorted[i] << '\t'
				<< outputData.velocity_zNonSorted[i] << '\n';

	output_velocityVector.close();

	for (std::size_t k = 0; k < fieldDataFileName_concentration.size(); ++k)
	{
		std::ofstream output_concentrationScalar {
				fieldDataFileName_concentration[k] };
		output_concentrationScalar.precision(20);
		if (!output_concentrationScalar.is_open())
			throw std::ofstream::failure(
					std::string("Couldn't create output file for field data ")
							+ std::string(fieldDataFileName_concentration[k])
							+ std::string("."));

		for (std::size_t i = 0; i < outputData.concentrationNonSorted[k].size();
				++i)
			output_concentrationScalar
					<< outputData.concentrationNonSorted[k][i] << '\n';

		output_concentrationScalar.close();
	}
}

void schemi::output::mixedZoneWidth1D(const structForOutput & outputData,
		const scalar Time)
{
	std::string timeWidthFileName { "./result/timeWidth.tsv" };

	std::ofstream timeWidthFile { timeWidthFileName, std::ios::app };
	timeWidthFile.precision(ioPrecision);

	if (timeWidthFile.is_open())
		std::cout << "./result/timeWidth.tsv is opened." << std::endl;
	else
		[[unlikely]]
		throw std::ofstream::failure("Couldn't create ./result/timeWidth.tsv.");

	scalar rL { 0 }, rR { 0 };
	bool isFounded { false };
	for (std::size_t i = 0; i < outputData.concentration[0].size() - 1; ++i)
	{
		const bool b1 { (outputData.concentration[1][i]
				/ outputData.concentration[0][i] >= zeroMix)
				&& (outputData.concentration[1][i + 1]
						/ outputData.concentration[0][i + 1] < zeroMix) };
		const bool b2 { (outputData.concentration[1][i]
				/ outputData.concentration[0][i] < zeroMix)
				&& (outputData.concentration[1][i + 1]
						/ outputData.concentration[0][i + 1] >= zeroMix) };

		if (b1 || b2)
		{
			rL = 0.5 * (outputData.x_coord[i] + outputData.x_coord[i + 1]);
			isFounded = true;
			break;
		}
	}

	if (!isFounded)
	{
		const std::size_t NCells { outputData.concentration[0].size() };
		if ((outputData.concentration[1][0] / outputData.concentration[0][0])
				< (outputData.concentration[1][NCells - 1]
						/ outputData.concentration[0][NCells - 1]))
			rL = outputData.x_coord[0];
		else
			rL = outputData.x_coord[NCells - 1];
	}

	isFounded = false;
	for (std::size_t i = 0; i < outputData.concentration[0].size() - 1; ++i)
	{
		const bool b1 { (outputData.concentration[1][i]
				/ outputData.concentration[0][i] > (1 - zeroMix))
				&& (outputData.concentration[1][i + 1]
						/ outputData.concentration[0][i + 1] <= (1 - zeroMix)) };
		const bool b2 { (outputData.concentration[1][i]
				/ outputData.concentration[0][i] <= (1 - zeroMix))
				&& (outputData.concentration[1][i + 1]
						/ outputData.concentration[0][i + 1] > (1 - zeroMix)) };

		if (b1 || b2)
		{
			rR = 0.5 * (outputData.x_coord[i] + outputData.x_coord[i + 1]);
			isFounded = true;
			break;
		}

	}

	if (!isFounded)
	{
		const std::size_t NCells { outputData.concentration[0].size() };
		if ((outputData.concentration[1][0] / outputData.concentration[0][0])
				> (outputData.concentration[1][NCells - 1]
						/ outputData.concentration[0][NCells - 1]))
			rR = outputData.x_coord[0];
		else
			rR = outputData.x_coord[NCells - 1];
	}

	timeWidthFile << Time << '\t' << std::abs(rR - rL) << '\n';

	timeWidthFile.close();
}
