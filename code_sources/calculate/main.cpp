#include "grid_field.h"
#include "sim_params.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <filesystem>


int main() {
    SimParams::writeParameters();

    GridField<SimParams::boundaryType> phi(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> eta(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> tempField(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> energies((SimParams::timeSteps / SimParams::outputInterval) + 1, 0.0, 0.0);

    // Установка начальных условий
    double waveNumber = 2 * M_PI * 4 / SimParams::domainLength;
    // phi.setSinusoidalInitialCondition(waveNumber, 0.4, SimParams::gridSpacing, M_PI/2);
    // phi.setBookInitialCondition(SimParams::domainLength, SimParams::gridSpacing);
    phi.setRandomInitialCondition(-0.2, 0.2);
    // phi.loadInitialConditionWithOffset(0.45, SimParams::Paths::initialConditionFile);

    double hSquared = SimParams::gridSpacing * SimParams::gridSpacing;
    // Основной цикл по времени
    for (int step = 0; step <= SimParams::timeSteps; ++step) {

        // Первый промежуточный шаг - вычисление eta
        for (int i = 0; i < SimParams::gridSize; ++i) {
            eta[i] = phi[i] + phi.laplacian(i, hSquared);
        }

        // Второй промежуточный шаг - вычисление tempField
        for (int i = 0; i < SimParams::gridSize; ++i) {
            tempField[i] = (-SimParams::epsilon * phi[i] + eta[i]
                            + eta.laplacian(i, hSquared)
                            + phi[i] * phi[i] * phi[i]);
        }

        // Финальный шаг - вычисление следующего значения phi
        for (int i = 0; i < SimParams::gridSize; ++i) {
            phiNext[i] = phi[i] + SimParams::timeStep * tempField.laplacian(i, hSquared);
        }

        // Обработка и сохранение результатов с заданной периодичностью
        if (step % SimParams::outputInterval == 0) {
            // Проверка значений на корректность
            phi.validateValues(10.0);

            // Вычисление энергии для текущего шага
            double energy = 0.0;
            for (int i = 0; i < SimParams::gridSize; ++i) {
                double term = tempField[i] * phi[i] / 2.0 - phi[i] * phi[i] * phi[i] * phi[i] / 4.0;
                energy += term * SimParams::gridSpacing;
            }
            energies[step / SimParams::outputInterval] = energy;

            // Сохранение текущего состояния
            phi.saveToFile(SimParams::Paths::getDataFilePath(step));

            // Отображение прогресса
            double progressPercent = 100.0 * step / SimParams::timeSteps;
            std::cout << progressPercent << "% - saved: " << SimParams::Paths::getDataFilePath(step) << std::endl;
        }

        // Обмен текущего и следующего полей phi
        swap(phi, phiNext);
    }

    // Сохранение итоговых данных энергии
    std::string energyFilePath = SimParams::Paths::getEnergiesFile();
    energies.saveToFile(energyFilePath);
    std::cout << "Energy data saved to: " << energyFilePath << std::endl;

    return 0;
}
