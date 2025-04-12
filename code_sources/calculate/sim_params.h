#ifndef SIM_PARAMS_H
#define SIM_PARAMS_H

#include <string>


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

    // File paths
    namespace Paths {
        // Output directory
        const std::string outputDir = "bin_data/";

        // Parameter file
        const std::string paramsFile = outputDir + "params.txt";

        // Data files
        const std::string getDataFilePath(int step) {
            return outputDir + "data_" + std::to_string(step) + ".bin";
        }

        // Energy data file
        const std::string energiesFile = outputDir + "energies.bin";

        // Initial condition file
        const std::string initialConditionFile = outputDir + "for_initial.bin";
    }
}

#endif // SIM_PARAMS_H
