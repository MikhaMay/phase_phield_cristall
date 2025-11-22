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
        mu[i] = phi[i]
              + 2.0 * phiLaplacian[i]
              + phiLaplacian.laplacian(i, hSquared)
              + phi[i] * phi[i] * phi[i]
              + r * phi[i];
}


// φ_{k+1/2} = (φ_k + φ_{k+1}) / 2
template <BoundaryType BType>
void calculateHalfAverages(
    const GridField<BType>& phi,
    GridField<BType>& phiHalf
) {
    for (int k = 0; k < phi.size(); ++k) {
        phiHalf[k] = 0.5 * (phi[k] + phi[k+1]);
    }
}

// (∇φ)_k = (φ_{k+1/2} - φ_{k-1/2}) / h
template <BoundaryType BType>
void calculateGradientFromHalfAtNodes(
    const GridField<BType>& phiHalf,
    GridField<BType>& gradPhi,
    double h
) {
    for (int k = 0; k < phiHalf.size(); ++k) {
        gradPhi[k] = (phiHalf[k] - phiHalf[k-1]) / h;
    }
}

// (Δφ)_{k+1/2} = (φ_{k+3/2} - 2φ_{k+1/2} + φ_{k-1/2}) / h^2
// Реализуем как лапласиан на "полуцелой" сетке
template <BoundaryType BType>
void calculateLaplacianOnHalf(
    const GridField<BType>& phiHalf,
    GridField<BType>& deltaPhiHalf,
    double hSquared
) {
    for (int k = 0; k < phiHalf.size(); ++k) {
        deltaPhiHalf[k] = (phiHalf[k+1] - 2.0 * phiHalf[k] + phiHalf[k-1]) / hSquared;
    }
}

// (∇Δφ)_k = ((Δφ)_{k+1/2} - (Δφ)_{k-1/2}) / h
template <BoundaryType BType>
void calculateGradLaplacianFromHalfAtNodes(
    const GridField<BType>& deltaPhiHalf,
    GridField<BType>& gradDeltaPhi,
    double h
) {
    for (int k = 0; k < deltaPhiHalf.size(); ++k) {
        gradDeltaPhi[k] = (deltaPhiHalf[k] - deltaPhiHalf[k-1]) / h;
    }
}

// ξ = μ ∇φ - [ (r+2)φ ∇φ - 2(Δφ) ∇φ + (Δφ) (∇Δφ) + (φ² - 1)φ ∇φ ]
template <BoundaryType BType>
void calculateXiAtNodes(
    const GridField<BType>& mu,              // μ at nodes
    const GridField<BType>& phi,             // φ at nodes
    const GridField<BType>& lapPhiNode,      // Δφ at nodes
    const GridField<BType>& gradPhiNode,     // ∇φ at nodes (через полуцелые)
    const GridField<BType>& gradDeltaPhiNode,// ∇Δφ at nodes (через полуцелые)
    GridField<BType>& xi,                    // result
    double r
) {
    for (int k = 0; k < phi.size(); ++k) {
        double phi_k = phi[k];
        double mu_k = mu[k];
        double grad_phi_k = gradPhiNode[k];
        double lap_phi_k  = lapPhiNode[k];
        double grad_lap_phi_k = gradDeltaPhiNode[k];

        double bracket =
            (r + 2.0) * phi_k * grad_phi_k
            - 2.0 * lap_phi_k * grad_phi_k
            + lap_phi_k * grad_lap_phi_k
            + (phi_k * phi_k - 1.0) * phi_k * grad_phi_k;

        xi[k] = mu_k * grad_phi_k - bracket;
    }
}

// Гауссовское сглаживание (coarse graining)
template <BoundaryType BType>
void applyCoarseGraining(
    const GridField<BType>& inField,
    GridField<BType>& outField,
    double h,
    double a0
) {
    const double norm_factor = 1.0 / (std::sqrt(2.0 * M_PI) * a0);
    const int N = inField.size();
    
    // Определяем радиус окна сглаживания (например, 3*a0)
    const int radius = static_cast<int>(std::ceil(3.0 * a0 / h));

    for (int k = 0; k < N; ++k) {
        double sum = 0.0;
        
        // Симметричное окно: от -radius до +radius
        for (int offset = -radius; offset <= radius; ++offset) {
            // Периодический индекс
            int j = (k + offset + N) % N;
            if (j < 0) j += N;  // дополнительная защита
            
            double dist = offset * h;
            double gauss = norm_factor * std::exp(-(dist * dist) / (2.0 * a0 * a0));
            sum += inField[j] * gauss * h;
        }
        outField[k] = sum;
    }
}

int main() {
    SimParams::writeParameters();

    // Поля на узлах
    GridField<SimParams::boundaryType> phi(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiLaplacian(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiLaplacianNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> mu(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> muNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> v(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> vNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> vLaplacian(SimParams::gridSize, 0.0, 0.0);

    // Поля на полуцелых индексах (размер такой же, индексы понимаем как k+1/2)
    GridField<SimParams::boundaryType> phiHalf(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> deltaPhiHalf(SimParams::gridSize, 0.0, 0.0);

    // Производные в узлах, построенные через полуцелые
    GridField<SimParams::boundaryType> gradPhiNode(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> gradDeltaPhiNode(SimParams::gridSize, 0.0, 0.0);

    // ξ и его сглаженная версия
    GridField<SimParams::boundaryType> xi(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> coarseXi(SimParams::gridSize, 0.0, 0.0);

    // Энергия
    GridField<SimParams::boundaryType> energies((SimParams::timeSteps / SimParams::outputInterval) + 1, 0.0, 0.0);

    // Начальные условия
    phi.loadInitialCondition(SimParams::Paths::initialConditionFile);
    v.setRandomInitialCondition(0.0, 0.0); // нулевая скорость

    const double h = SimParams::gridSpacing;
    const double hSquared = h * h;
    const double dt = SimParams::timeStep;

    for (int step = 0; step <= SimParams::timeSteps; ++step) {
        // Текущие Δφ и μ
        calculateLaplacian(phi, phiLaplacian, hSquared);
        calculateChemicalPotential(phi, phiLaplacian, mu, hSquared, SimParams::r);

        // Шаг по φ (явный, перенос upwind)
        for (int i = 0; i < SimParams::gridSize; ++i) {
            double convection_term = 0.0;
            if (v[i] > 0.0)
                convection_term = v[i] * phi.nabla(i, h, false); // (φ_i - φ_{i-1})/h
            else
                convection_term = v[i] * phi.nabla(i, h, true);  // (φ_{i+1} - φ_i)/h

            phiNext[i] = phi[i] + dt * (SimParams::Gamma * mu.laplacian(i, hSquared) - convection_term);
        }

        // Поля для n+1
        calculateLaplacian(phiNext, phiLaplacianNext, hSquared);
        calculateChemicalPotential(phiNext, phiLaplacianNext, muNext, hSquared, SimParams::r);

        // Полуцелые: φ_{k+1/2}, затем ∇φ|_k
        calculateHalfAverages(phiNext, phiHalf);
        calculateGradientFromHalfAtNodes(phiHalf, gradPhiNode, h);

        // Δφ на полуцелых и затем ∇Δφ|_k
        calculateLaplacianOnHalf(phiHalf, deltaPhiHalf, hSquared);
        calculateGradLaplacianFromHalfAtNodes(deltaPhiHalf, gradDeltaPhiNode, h);

        // ξ^{n+1} в узлах
        calculateXiAtNodes(muNext, phiNext, phiLaplacianNext, gradPhiNode, gradDeltaPhiNode, xi, SimParams::r);

        // Гауссово усреднение ξ^{n+1}
        applyCoarseGraining(xi, coarseXi, h, SimParams::a_0);
        calculateLaplacian(v, vLaplacian, hSquared);

        // Обновление v: v^{n+1} = v^n + dt/rho0 * ( <ξ^{n+1}> + Γ_S Δ v^n )
        for (int i = 0; i < SimParams::gridSize; ++i) {
            vNext[i] = v[i] + dt * (coarseXi[i] + SimParams::Gamma_S * vLaplacian[i]) / SimParams::rho_0;
        }

        // Вывод/энергия
        if (step % SimParams::outputInterval == 0) {
            phi.validateValues(10.0);

            // Для энергии используем согласованный ∇φ через полуцелые для текущего φ
            calculateHalfAverages(phi, phiHalf);
            calculateGradientFromHalfAtNodes(phiHalf, gradPhiNode, h);

            double energy = 0.0;
            for (int i = 0; i < SimParams::gridSize; ++i) {
                double f = (SimParams::r + 2.0) * phi[i] * phi[i] / 2.0
                         - gradPhiNode[i] * gradPhiNode[i]
                         + phiLaplacian[i] * phiLaplacian[i] / 2.0
                         + (phi[i] * phi[i] - 1.0) * (phi[i] * phi[i] - 1.0) / 4.0;
                energy += f * h;
            }
            energies[step / SimParams::outputInterval] = energy;

            // Сохранение: tmp теперь — это xi (для ясности можно переименовать путь)
            phi.saveToFile(SimParams::Paths::getDataFilePath(step));
            v.saveToFile(SimParams::Paths::getVelocityFilePath(step));
            vLaplacian.saveToFile(SimParams::Paths::getVLaplacianFilePath(step));
            xi.saveToFile(SimParams::Paths::getXiFilePath(step));
            coarseXi.saveToFile(SimParams::Paths::getCoarseXiFilePath(step));
            mu.saveToFile(SimParams::Paths::getMuFilePath(step));

            double progressPercent = 100.0 * step / SimParams::timeSteps;
            std::cout << progressPercent << "% - saved: " << SimParams::Paths::getDataFilePath(step) << std::endl;
        }

        // swap
        swap(phi, phiNext);
        swap(v, vNext);
    }

    // Сохранение энергий
    std::string energyFilePath = SimParams::Paths::getEnergiesFile();
    energies.saveToFile(energyFilePath);
    std::cout << energies.size() << " energy values" << std::endl;
    std::cout << "Energy data saved to: " << energyFilePath << std::endl;

    return 0;
}
