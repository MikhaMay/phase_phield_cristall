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
    constexpr int gridSize = 100;
    constexpr double gridSpacing = domainLength / gridSize;
    constexpr double timeStep = 1e-7;
    constexpr int timeSteps = 5'000'000;
    // constexpr double epsilon = 0.4;
    constexpr int outputInterval = 50'000;
    constexpr double totalTime = timeStep * timeSteps;
    constexpr BoundaryType boundaryType = BoundaryType::Periodic;

    // Model parameters
    constexpr double r = -0.4;  // параметр из уравнения свободной энергии
    constexpr double Gamma = 1.0;  // подвижность
    constexpr double Gamma_S = 0.001;  // параметр диссипации скорости
    constexpr double rho_0 = 1; //0.001;  // плотность
    constexpr double a_0 = 80.0 / 14.0;  // постоянная кристаллической решетки

    // Custom tag for the current simulation
    const std::string simulationTag = "sHPFC_0_boundary_velocity";

    // File paths
    namespace Paths {
        // Base directories
        const std::string outputDir = "calculate_output_data/";
        const std::string inputDir = "calculate_input_data/";

        // Generate descriptive filename prefix with parameters
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

        // Parameter file
        const std::string paramsFile = outputDir + "params.yaml";

        // Initial condition file
        const std::string initialConditionFile = inputDir + "initial_condition.bin";

        // Data files with descriptive names
        inline std::string getDataFilePath(int step) {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "/" + std::to_string(step) + ".bin";
        }

        // Velocity data files with descriptive names
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

        // Energy data file with descriptive name
        inline std::string getEnergiesFile() {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "/energies.bin";
        }
    }

    void writeParameters() {
        std::string directoryName = Paths::outputDir + Paths::getDescriptivePrefix();
        std::filesystem::create_directories(directoryName);

        // Write parameters to the new directory
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
