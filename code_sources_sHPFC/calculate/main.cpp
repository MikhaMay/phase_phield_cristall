#include "grid_field.h"
#include "sim_params.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <cmath>


template <BoundaryType BType>
void calculateLaplacian(
    const GridField<BType>& phi,
    GridField<BType>& phiLaplacian,
    double hSquared
) {
    for (int i = 0; i < phi.size(); ++i)
        phiLaplacian[i] = phi.laplacian(i, hSquared);
}


template <BoundaryType BType>
void calculateChemicalPotential(
    const GridField<BType>& phi,
    const GridField<BType>& phiLaplacian,
    GridField<BType>& mu,
    double hSquared,
    double r
) {
    // μ = (1 + Δ)²φ + φ³ + rφ = φ + 2Δφ + Δ²φ + φ³ + rφ
    for (int i = 0; i < phi.size(); ++i)
        mu[i] = phi[i] + 2 * phiLaplacian[i] + phiLaplacian.laplacian(i, hSquared) + phi[i] * phi[i] * phi[i] + r * phi[i];
}


template <BoundaryType BType>
void calculateAB(
    const GridField<BType>& phiLaplacian,
    GridField<BType>& A,
    GridField<BType>& B,
    double hSquared
) {
    for (int i = 0; i < phiLaplacian.size(); ++i) {
        // B = Δφ
        B[i] = phiLaplacian[i];

        // A = 4Δφ + Δ²φ
        A[i] = 4 * phiLaplacian[i] + phiLaplacian.laplacian(i, hSquared);
    }
}

// Локальное tmp
template <BoundaryType BType>
void calculateTmp(
    const GridField<BType>& A,
    const GridField<BType>& B,
    const GridField<BType>& phi,
    GridField<BType>& tmp,
    double h
) {
    for (int i = 0; i < phi.size(); ++i) {
        // A * nabla(phi)
        double A_term = 0.0;
        if (A[i] > 0)
            A_term = A[i] * phi.nabla(i, h, false);
        else
            A_term = A[i] * phi.nabla(i, h, true);

        // B * nabla(B)
        double B_term = 0.0;
        if (B[i] > 0)
            B_term = B[i] * B.nabla(i, h, false);
        else
            B_term = B[i] * B.nabla(i, h, true);

        // Итоговое значение tmp
        tmp[i] = A_term - B_term;
    }
}

// Гауссовское сглаживание (coarse graining)
template <BoundaryType BType>
void applyCoarseGraining(
    const GridField<BType>& tmp,
    GridField<BType>& coarseTmp,
    double h,
    double a0
) {
    const double norm_factor = 1.0 / (std::sqrt(2.0 * M_PI * a0 * a0));

    for (int k = 0; k < tmp.size(); ++k) {
        double sum = 0.0;

        // Применяем гауссово сглаживание по всем точкам
        for (int j = 0; j < tmp.size(); ++j) {
            double dist = (k - j) * h;
            double gauss = norm_factor * std::exp(-(dist * dist) / (2.0 * a0 * a0));
            sum += tmp[j] * gauss * h;
        }

        coarseTmp[k] = sum;
    }
}

int main() {
    SimParams::writeParameters();

    // Создание сеточных полей
    GridField<SimParams::boundaryType> phi(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiLaplacian(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiLaplacianNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> mu(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> v(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> vNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> A(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> B(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> tmp(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> coarseTmp(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> energies((SimParams::timeSteps / SimParams::outputInterval) + 1, 0.0, 0.0);

    // Установка начальных условий
    // phi.setRandomInitialCondition(-0.2, 0.2);
    // v.setRandomInitialCondition(-0.01, 0.01);  // Начальное поле скоростей

    phi.loadInitialCondition(SimParams::Paths::initialConditionFile);
    v.setRandomInitialCondition(0.0, 0.0);

    double h = SimParams::gridSpacing;
    double hSquared = h * h;
    double dt = SimParams::timeStep;

    // Основной цикл по времени
    for (int step = 0; step <= SimParams::timeSteps; ++step) {
        calculateLaplacian(phi, phiLaplacian, hSquared);
        // Вычисление химического потенциала
        calculateChemicalPotential(phi, phiLaplacian, mu, hSquared, SimParams::r);

        // Вычисление phi на следующем временном шаге
        for (int i = 0; i < SimParams::gridSize; ++i) {
            // Конвективный член с upwind схемой
            double convection_term = 0.0;
            if (v[i] > 0)
                convection_term = v[i] * phi.nabla(i, h, false);
            else
                convection_term = v[i] * phi.nabla(i, h, true);

            // Обновление phi
            phiNext[i] = phi[i] + dt * (SimParams::Gamma * mu.laplacian(i, hSquared) - convection_term);
        }


        calculateLaplacian(phiNext, phiLaplacianNext, hSquared);
        // Расчет параметров A и B
        calculateAB(phiLaplacianNext, A, B, hSquared);
        // Расчет локального tmp
        calculateTmp(A, B, phiNext, tmp, h);
        // Применение гауссовского сглаживания
        applyCoarseGraining(tmp, coarseTmp, h, SimParams::a_0);
        // Обновление скорости
        for (int i = 0; i < SimParams::gridSize; ++i) {
            vNext[i] = v[i] + dt * (coarseTmp[i] + SimParams::Gamma_S * v.laplacian(i, hSquared)) / SimParams::rho_0;
        }

        // Обработка и сохранение результатов с заданной периодичностью
        if (step % SimParams::outputInterval == 0) {
            // Проверка значений на корректность
            phi.validateValues(10.0);
            // v.validateValues(10.0);

            // Вычисление энергии для текущего шага
            double energy = 0.0;
            for (int i = 0; i < SimParams::gridSize; ++i) {
                // Свободная энергия: f = (r+2)φ²/2 - (∇φ)²/2 + (Δφ)²/2 + (φ²-1)²/4
                double grad_phi = phi.nabla(i, h);

                double f = (SimParams::r + 2) * phi[i] * phi[i] / 2.0 
                         - grad_phi * grad_phi / 2.0 
                         + phiLaplacian[i] * phiLaplacian[i] / 2.0 
                         + (phi[i] * phi[i] - 1.0) * (phi[i] * phi[i] - 1.0) / 4.0;

                energy += f * h;
            }
            energies[step / SimParams::outputInterval] = energy;

            // Сохранение текущего состояния
            phi.saveToFile(SimParams::Paths::getDataFilePath(step));
            v.saveToFile(SimParams::Paths::getVelocityFilePath(step));

            // Отображение прогресса
            double progressPercent = 100.0 * step / SimParams::timeSteps;
            std::cout << progressPercent << "% - saved: " << SimParams::Paths::getDataFilePath(step) << std::endl;
        }

        // Обмен текущего и следующего полей
        swap(phi, phiNext);
        swap(v, vNext);
    }

    // Сохранение итоговых данных энергии
    std::string energyFilePath = SimParams::Paths::getEnergiesFile();
    energies.saveToFile(energyFilePath);
    std::cout << "Energy data saved to: " << energyFilePath << std::endl;

    return 0;
}
