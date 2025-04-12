#include "grid_field.h"
#include "sim_params.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <filesystem>


void writeParameters(const std::string& path) {
    std::ofstream os(path);
    if (!os) {
        throw std::runtime_error("Failed to open parameters file for writing: " + path);
    }

    os << SimParams::timeSteps << std::endl;
    os << SimParams::outputInterval << std::endl;
    os << SimParams::epsilon << std::endl;

    os.close();
    if (os.fail()) {
        throw std::runtime_error("Error closing file after writing parameters: " + path);
    }
}

int main() {
    writeParameters(SimParams::Paths::paramsFile);
    
    // Инициализация полей с периодическими граничными условиями
    GridField<BoundaryType::Periodic> phi(SimParams::gridSize, 0.0, 0.0);
    GridField<BoundaryType::Periodic> phiNext(SimParams::gridSize, 0.0, 0.0);
    GridField<BoundaryType::Periodic> eta(SimParams::gridSize, 0.0, 0.0);
    GridField<BoundaryType::Periodic> tempField(SimParams::gridSize, 0.0, 0.0);
    GridField<BoundaryType::Periodic> energies((SimParams::timeSteps / SimParams::outputInterval) + 1, 0.0, 0.0);

    // Установка начальных условий
    // phi.setSinusoidalInitialCondition(1.2, 0.4, SimParams::gridSpacing, SimParams::PI/2);
    // phi.setBookInitialCondition(SimParams::domainLength, SimParams::gridSpacing);
    phi.setRandomInitialCondition(-0.2, 0.2);
    // phi.loadInitialConditionWithOffset(0.45, SimParams::Paths::initialConditionFile);

    // Основной цикл по времени
    for (int step = 0; step <= SimParams::timeSteps; ++step) {
        double hSquared = SimParams::gridSpacing * SimParams::gridSpacing;
        
        // Первый промежуточный шаг - вычисление eta
        for (int i = 0; i < SimParams::gridSize; ++i) {
            eta[i] = phi[i] + (phi[i-1] - 2 * phi[i] + phi[i+1]) / hSquared;
        }

        // Второй промежуточный шаг - вычисление tempField
        for (int i = 0; i < SimParams::gridSize; ++i) {
            tempField[i] = (-SimParams::epsilon * phi[i] + eta[i] + 
                            (eta[i-1] - 2 * eta[i] + eta[i+1]) / hSquared + 
                            phi[i] * phi[i] * phi[i]);
        }

        // Финальный шаг - вычисление следующего значения phi
        for (int i = 0; i < SimParams::gridSize; ++i) {
            phiNext[i] = phi[i] + SimParams::timeStep * (tempField[i-1] - 2 * tempField[i] + tempField[i+1]) / hSquared;
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
            std::cout << progressPercent << "%" << std::endl;
        }

        // Обмен текущего и следующего полей phi
        swap(phi, phiNext);
    }

    // Сохранение итоговых данных энергии
    energies.saveToFile(SimParams::Paths::energiesFile);

    return 0;
}
