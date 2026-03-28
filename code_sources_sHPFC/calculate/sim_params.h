#ifndef SIM_PARAMS_H
#define SIM_PARAMS_H

#include "grid_field.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <stdexcept>

namespace SimParams {
    // Simulation parameters
    constexpr double domainLength = 80.0;
    constexpr int gridSize = 101; // теперь считаем, что узлы включают обе границы
    constexpr double gridSpacing = domainLength / (gridSize - 1);
    constexpr double timeStep = 1e-6;
    constexpr int timeSteps = 5'000'000;
    constexpr int outputInterval = 50'000;
    constexpr double totalTime = timeStep * timeSteps;
    constexpr BoundaryType boundaryType = BoundaryType::Fix;

    // Model parameters
    constexpr double r = -0.4;
    constexpr double Gamma = 1.0;
    constexpr double Gamma_S = 1.0;
    constexpr double rho_0 = 0.001;
    constexpr double a_0 = 80.0 / 14.0;

    const std::string simulationTag = "sHPFC_spring_in_jar_01";

    namespace Paths {
        const std::string outputDir = "calculate_output_data/";
        const std::string inputDir = "calculate_input_data/";

        inline std::string getDescriptivePrefix() {
            std::ostringstream prefix;
            prefix << "L" << std::fixed << std::setprecision(1) << domainLength
                   << "_N" << gridSize
                   << "_r" << std::setprecision(2) << r
                   << "_dt" << std::scientific << std::setprecision(1) << timeStep
                   << "_" << SimParams::boundaryType
                   << "_" << simulationTag;
            return prefix.str();
        }

        const std::string paramsFile = outputDir + "params.yaml";
        const std::string initialConditionFile = inputDir + "initial_condition.bin";

        inline std::string getDataFilePath(int step) {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "/" + std::to_string(step) + ".bin";
        }

        inline std::string getVelocityFilePath(int step) {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "/v_" + std::to_string(step) + ".bin";
        }

        inline std::string getVLaplacianFilePath(int step) {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "/vlap_" + std::to_string(step) + ".bin";
        }

        inline std::string getXiFilePath(int step) {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "/xi_" + std::to_string(step) + ".bin";
        }

        inline std::string getCoarseXiFilePath(int step) {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "/coarse_xi_" + std::to_string(step) + ".bin";
        }

        inline std::string getMuFilePath(int step) {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "/mu_" + std::to_string(step) + ".bin";
        }

        inline std::string getEnergiesFile() {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "/energies.bin";
        }
    }

    void writeParameters() {
        std::string directoryName = Paths::outputDir + Paths::getDescriptivePrefix();
        std::filesystem::create_directories(directoryName);

        std::string paramsFilePath = directoryName + "/params.yaml";
        std::ofstream paramsFile(paramsFilePath);
        if (!paramsFile) {
            throw std::runtime_error("Failed to open parameters file for writing: " + paramsFilePath);
        }

        paramsFile << "domainLength: " << SimParams::domainLength << std::endl;
        paramsFile << "gridSize: " << SimParams::gridSize << std::endl;
        paramsFile << "gridSpacing: " << SimParams::gridSpacing << std::endl;
        paramsFile << "timeStep: " << SimParams::timeStep << std::endl;
        paramsFile << "timeSteps: " << SimParams::timeSteps << std::endl;
        paramsFile << "outputInterval: " << SimParams::outputInterval << std::endl;
        paramsFile << "r: " << SimParams::r << std::endl;
        paramsFile << "Gamma: " << SimParams::Gamma << std::endl;
        paramsFile << "Gamma_S: " << SimParams::Gamma_S << std::endl;
        paramsFile << "rho_0: " << SimParams::rho_0 << std::endl;
        paramsFile << "a_0: " << SimParams::a_0 << std::endl;
        paramsFile << "boundaryType: " << SimParams::boundaryType << std::endl;
        paramsFile << "simulationTag: " << SimParams::simulationTag << std::endl;

        paramsFile.close();
        if (paramsFile.fail()) {
            throw std::runtime_error("Error closing file after writing parameters: " + paramsFilePath);
        }

        std::string latestRunFilePath = Paths::outputDir + "/latest_run.txt";
        std::ofstream latestRunFile(latestRunFilePath);
        if (!latestRunFile) {
            throw std::runtime_error("Failed to open latest run file for writing: " + latestRunFilePath);
        }

        latestRunFile << directoryName << std::endl;

        latestRunFile.close();
        if (latestRunFile.fail()) {
            throw std::runtime_error("Error closing latest run file: " + latestRunFilePath);
        }
    }
}

#endif // SIM_PARAMS_H
