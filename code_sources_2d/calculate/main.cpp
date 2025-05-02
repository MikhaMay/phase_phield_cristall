#include "grid_field.h"
#include "sim_params.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <filesystem>

int main() {
    SimParams::writeParameters();

    GridField<SimParams::boundaryType> phi(SimParams::gridSizeX, SimParams::gridSizeY, 0.0);
    GridField<SimParams::boundaryType> phiNext(SimParams::gridSizeX, SimParams::gridSizeY, 0.0);
    GridField<SimParams::boundaryType> eta(SimParams::gridSizeX, SimParams::gridSizeY, 0.0);
    GridField<SimParams::boundaryType> tempField(SimParams::gridSizeX, SimParams::gridSizeY, 0.0);
    GridField<SimParams::boundaryType> energies((SimParams::timeSteps / SimParams::outputInterval) + 1, 1, 0.0);

    // Set initial conditions
    // double waveNumberX = 2 * M_PI * 4 / SimParams::domainLengthX;
    // double waveNumberY = 2 * M_PI * 4 / SimParams::domainLengthY;

    // Choose one of the following initial conditions:
    phi.setLoadHalfFromFileInitialCondition(SimParams::Paths::initialConditionFile, -0.9);
    // phi.setSinusoidalInitialCondition(waveNumberX, waveNumberY, 0.4, SimParams::gridSpacingX, SimParams::gridSpacingY, M_PI/2, M_PI/2);
    // phi.setBookInitialCondition(SimParams::domainLengthX, SimParams::domainLengthY, SimParams::gridSpacingX, SimParams::gridSpacingY);
    // phi.loadInitialConditionWithOffset(0.45, SimParams::Paths::initialConditionFile);

    // Main time loop
    for (unsigned long step = 0; step <= SimParams::timeSteps; ++step) {
        double hxSquared = SimParams::gridSpacingX * SimParams::gridSpacingX;
        double hySquared = SimParams::gridSpacingY * SimParams::gridSpacingY;

        // First intermediate step - calculate eta
        for (int i = 0; i < SimParams::gridSizeX; ++i) {
            for (int j = 0; j < SimParams::gridSizeY; ++j) {
                eta(i, j) = phi(i, j) + phi.laplacian(i, j, hxSquared, hySquared);
            }
        }

        // Second intermediate step - calculate tempField
        for (int i = 0; i < SimParams::gridSizeX; ++i) {
            for (int j = 0; j < SimParams::gridSizeY; ++j) {
                tempField(i, j) = (-SimParams::epsilon * phi(i, j) + eta(i, j)
                              + eta.laplacian(i, j, hxSquared, hySquared)
                              + phi(i, j) * phi(i, j) * phi(i, j));
            }
        }

        // Final step - calculate next value of phi
        for (int i = 0; i < SimParams::gridSizeX; ++i) {
            for (int j = 0; j < SimParams::gridSizeY; ++j) {
                phiNext(i, j) = phi(i, j) + SimParams::timeStep * tempField.laplacian(i, j, hxSquared, hySquared);
            }
        }

        // Process and save results at specified intervals
        if (step % SimParams::outputInterval == 0) {
            // Validate values
            phi.validateValues(10.0);

            // Calculate energy for the current step
            double energy = 0.0;
            for (int i = 0; i < SimParams::gridSizeX; ++i) {
                for (int j = 0; j < SimParams::gridSizeY; ++j) {
                    double term = tempField(i, j) * phi(i, j) / 2.0 - phi(i, j) * phi(i, j) * phi(i, j) * phi(i, j) / 4.0;
                    energy += term * SimParams::gridSpacingX * SimParams::gridSpacingY;
                }
            }
            energies(step / SimParams::outputInterval, 0) = energy;

            // Save current state
            phi.saveToFile(SimParams::Paths::getDataFilePath(step));

            // Show progress
            double progressPercent = 100.0 * step / SimParams::timeSteps;
            std::cout << progressPercent << "% - saved: " << SimParams::Paths::getDataFilePath(step) << std::endl;
        }

        // Swap current and next phi fields
        swap(phi, phiNext);
    }

    // Save final energy data
    std::string energyFilePath = SimParams::Paths::getEnergiesFile();
    energies.saveToFile(energyFilePath);
    std::cout << "Energy data saved to: " << energyFilePath << std::endl;

    return 0;
}
