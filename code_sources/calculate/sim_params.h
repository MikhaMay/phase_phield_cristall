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
    constexpr int gridSize = 200;
    constexpr double gridSpacing = domainLength / gridSize;
    constexpr double timeStep = 1e-4;
    constexpr int timeSteps = 1'000'000;
    constexpr double epsilon = 0.4;
    constexpr int outputInterval = 10'000;
    constexpr double totalTime = timeStep * timeSteps;
    constexpr double PI = 3.14159265358979323846;
    // Тип граничных условий для симуляции
    constexpr BoundaryType boundaryType = BoundaryType::Periodic;

    // Custom tag for the current simulation
    const std::string simulationTag = "default";

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
                   << "_eps" << std::setprecision(2) << epsilon
                   << "_dt" << std::scientific << std::setprecision(1) << timeStep
                   << "_" << simulationTag;
            return prefix.str();
        }

        // Parameter file
        const std::string paramsFile = outputDir + "params.yaml";

        // Initial condition file - now in input directory
        const std::string initialConditionFile = inputDir + "initial_condition.bin";

        // Custom initial condition paths
        inline std::string getInitialConditionPath(const std::string& name) {
            return inputDir + name + ".bin";
        }

        // Data files with descriptive names
        inline std::string getDataFilePath(int step) {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "_step" + std::to_string(step) + ".bin";
        }

        // Energy data file with descriptive name
        inline std::string getEnergiesFile() {
            std::string prefix = getDescriptivePrefix();
            return outputDir + prefix + "_energies.bin";
        }
    }

    void writeParameters(const std::string& path) {
        std::ofstream os(path);
        if (!os) {
            throw std::runtime_error("Failed to open parameters file for writing: " + path);
        }

        os << "domainLength: " << SimParams::domainLength << std::endl;
        os << "gridSize: " << SimParams::gridSize << std::endl;
        os << "gridSpacing: " << SimParams::gridSpacing << std::endl;
        os << "timeStep: " << SimParams::timeStep << std::endl;
        os << "timeSteps: " << SimParams::timeSteps << std::endl;
        os << "outputInterval: " << SimParams::outputInterval << std::endl;
        os << "epsilon: " << SimParams::epsilon << std::endl;
        os << "boundaryType: " << (SimParams::boundaryType == BoundaryType::Periodic ? "Periodic" : "Fixed") << std::endl;
        os << "simulationTag: " << SimParams::simulationTag << std::endl;
    
        os.close();
        if (os.fail()) {
            throw std::runtime_error("Error closing file after writing parameters: " + path);
        }
    }
}

#endif // SIM_PARAMS_H
