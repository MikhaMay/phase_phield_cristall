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
    // μ = (1 + Δ)^2 φ + φ^3 + rφ = φ + 2Δφ + Δ^2φ + φ^3 + rφ
    for (int i = 0; i < phi.size(); ++i)
        mu[i] = phi[i]
              + 2.0 * phiLaplacian[i]
              + phiLaplacian.laplacian(i, hSquared)
              + phi[i] * phi[i] * phi[i]
              + r * phi[i];
}

template <BoundaryType BType>
void calculateHalfAverages(
    const GridField<BType>& phi,
    GridField<BType>& phiHalf
) {
    for (int k = 0; k < phi.size(); ++k) {
        phiHalf[k] = 0.5 * (phi[k] + phi[k + 1]);
    }
}

template <BoundaryType BType>
void calculateGradientFromHalfAtNodes(
    const GridField<BType>& phiHalf,
    GridField<BType>& gradPhi,
    double h
) {
    for (int k = 0; k < phiHalf.size(); ++k) {
        gradPhi[k] = (phiHalf[k] - phiHalf[k - 1]) / h;
    }
}

template <BoundaryType BType>
void calculateLaplacianOnHalf(
    const GridField<BType>& phiHalf,
    GridField<BType>& deltaPhiHalf,
    double hSquared
) {
    for (int k = 0; k < phiHalf.size(); ++k) {
        deltaPhiHalf[k] = (phiHalf[k + 1] - 2.0 * phiHalf[k] + phiHalf[k - 1]) / hSquared;
    }
}

template <BoundaryType BType>
void calculateGradLaplacianFromHalfAtNodes(
    const GridField<BType>& deltaPhiHalf,
    GridField<BType>& gradDeltaPhi,
    double h
) {
    for (int k = 0; k < deltaPhiHalf.size(); ++k) {
        gradDeltaPhi[k] = (deltaPhiHalf[k] - deltaPhiHalf[k - 1]) / h;
    }
}

template <BoundaryType BType>
void calculateXiAtNodes(
    const GridField<BType>& mu,
    const GridField<BType>& phi,
    const GridField<BType>& lapPhiNode,
    const GridField<BType>& gradPhiNode,
    const GridField<BType>& gradDeltaPhiNode,
    GridField<BType>& xi,
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

template <BoundaryType BType>
void applyCoarseGraining(
    const GridField<BType>& inField,
    GridField<BType>& outField,
    double h,
    double a0
) {
    const double norm_factor = 1.0 / (std::sqrt(2.0 * M_PI) * a0);
    const int N = inField.size();
    const int radius = static_cast<int>(std::ceil(3.0 * a0 / h));

    for (int k = 0; k < N; ++k) {
        double sum = 0.0;

        for (int offset = -radius; offset <= radius; ++offset) {
            int j = (k + offset + N) % N;
            if (j < 0) j += N;

            double dist = offset * h;
            double gauss = norm_factor * std::exp(-(dist * dist) / (2.0 * a0 * a0));
            sum += inField[j] * gauss * h;
        }
        outField[k] = sum;
    }
}

template <BoundaryType BType>
void setZeroInitialCondition(GridField<BType>& field) {
    for (int i = 0; i < field.size(); ++i) {
        field[i] = 0.0;
    }
}

// f_ext(x,t) = A sin(2*pi*m*x/L) cos(omega*t)
inline double externalForce(double x, double t) {
    if (!SimParams::useExternalForce) {
        return 0.0;
    }

    return SimParams::forceAmplitude;

    // const double q = 2.0 * M_PI * SimParams::forceMode / SimParams::domainLength;
    // return SimParams::forceAmplitude * std::sin(q * x) * std::cos(SimParams::forceOmega * t);
}

template <BoundaryType BType>
void calculateExternalForceField(
    GridField<BType>& fExtField,
    double currentTime,
    int step,
    double h
) {
    setZeroInitialCondition(fExtField);

    if (!SimParams::useExternalForce) {
        return;
    }

    if (step >= 50000) {
        return;
    }

    const double limit = step / 5000.0;

    for (int i = 0; i < fExtField.size(); ++i) {
        const double x = i * h;
        const double force = externalForce(x, currentTime);

        if (i < limit) {
            fExtField[i] = force;
        }
        else if (i > SimParams::gridSize - 1 - limit) {
            fExtField[i] = -force;
        }
    }
}

int main() {
    SimParams::writeParameters();

    GridField<SimParams::boundaryType> phi(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiLaplacian(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> phiLaplacianNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> mu(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> muNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> v(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> vNext(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> vLaplacian(SimParams::gridSize, 0.0, 0.0);

    GridField<SimParams::boundaryType> phiHalf(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> deltaPhiHalf(SimParams::gridSize, 0.0, 0.0);

    GridField<SimParams::boundaryType> gradPhiNode(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> gradDeltaPhiNode(SimParams::gridSize, 0.0, 0.0);

    GridField<SimParams::boundaryType> xi(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> coarseXi(SimParams::gridSize, 0.0, 0.0);
    GridField<SimParams::boundaryType> fExtField(SimParams::gridSize, 0.0, 0.0);

    GridField<SimParams::boundaryType> energies(
        (SimParams::timeSteps / SimParams::outputInterval) + 1, 0.0, 0.0
    );

    phi.loadInitialCondition(SimParams::Paths::initialConditionFile);
    setZeroInitialCondition(v);

    const double h = SimParams::gridSpacing;
    const double hSquared = h * h;
    const double dt = SimParams::timeStep;

    for (int step = 0; step <= SimParams::timeSteps; ++step) {
        const double currentTime = step * dt;

        calculateLaplacian(phi, phiLaplacian, hSquared);
        calculateChemicalPotential(phi, phiLaplacian, mu, hSquared, SimParams::r);

        // Шаг по φ: PFC + перенос Lax–Wendroff
        for (int i = 0; i < SimParams::gridSize; ++i) {
            double v_loc = v[i];
            double lambda = v_loc * dt / h;

            double conv =
                -0.5 * lambda * (phi[i + 1] - phi[i - 1]) +
                 0.5 * lambda * lambda * (phi[i + 1] - 2.0 * phi[i] + phi[i - 1]);

            phiNext[i] = phi[i]
                       + dt * (SimParams::Gamma * mu.laplacian(i, hSquared))
                       + conv;
        }

        calculateLaplacian(phiNext, phiLaplacianNext, hSquared);
        calculateChemicalPotential(phiNext, phiLaplacianNext, muNext, hSquared, SimParams::r);

        calculateHalfAverages(phiNext, phiHalf);
        calculateGradientFromHalfAtNodes(phiHalf, gradPhiNode, h);

        calculateLaplacianOnHalf(phiHalf, deltaPhiHalf, hSquared);
        calculateGradLaplacianFromHalfAtNodes(deltaPhiHalf, gradDeltaPhiNode, h);

        calculateXiAtNodes(
            muNext,
            phiNext,
            phiLaplacianNext,
            gradPhiNode,
            gradDeltaPhiNode,
            xi,
            SimParams::r
        );

        applyCoarseGraining(xi, coarseXi, h, SimParams::a_0);
        calculateLaplacian(v, vLaplacian, hSquared);
        calculateExternalForceField(fExtField, currentTime, step, h);

        // Обновление v с внешней силой
        for (int i = 0; i < SimParams::gridSize; ++i) {
            vNext[i] = v[i]
                     + dt * (coarseXi[i] + SimParams::Gamma_S * vLaplacian[i] + fExtField[i])
                            / SimParams::rho_0;
        }

        if (step % SimParams::outputInterval == 0) {
            phi.validateValues(10.0);

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

            phi.saveToFile(SimParams::Paths::getDataFilePath(step));
            v.saveToFile(SimParams::Paths::getVelocityFilePath(step));
            vLaplacian.saveToFile(SimParams::Paths::getVLaplacianFilePath(step));
            xi.saveToFile(SimParams::Paths::getXiFilePath(step));
            coarseXi.saveToFile(SimParams::Paths::getCoarseXiFilePath(step));
            mu.saveToFile(SimParams::Paths::getMuFilePath(step));
            fExtField.saveToFile(SimParams::Paths::getFExtFilePath(step));

            double progressPercent = 100.0 * step / SimParams::timeSteps;
            std::cout << progressPercent << "% - saved: "
                      << SimParams::Paths::getDataFilePath(step) << std::endl;
        }

        swap(phi, phiNext);
        swap(v, vNext);
    }

    std::string energyFilePath = SimParams::Paths::getEnergiesFile();
    energies.saveToFile(energyFilePath);
    std::cout << energies.size() << " energy values" << std::endl;
    std::cout << "Energy data saved to: " << energyFilePath << std::endl;

    return 0;
}
