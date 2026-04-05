#ifndef GRID_FIELD_CPP
#define GRID_FIELD_CPP

#include "grid_field.h"

std::ostream& operator<<(std::ostream& os, const BoundaryType& type) {
    switch (type) {
        case BoundaryType::Periodic:
            os << "Periodic";
            break;
        case BoundaryType::Fix:
            os << "Fix";
            break;
        default:
            throw std::invalid_argument("Unknown BoundaryType value");
    }
    return os;
}

template <BoundaryType BType>
void GridField<BType>::saveToFile(const std::string& path) {
    // Open file for binary writing
    std::ofstream outfile(path, std::ios::binary);
    if (!outfile) {
        throw std::runtime_error("Failed to open file for writing: " + path);
    }

    // Write the size of the array
    outfile.write(reinterpret_cast<const char*>(&size_), sizeof(size_));

    // Write the array data
    outfile.write(reinterpret_cast<const char*>(values_), size_ * sizeof(double));
    
    outfile.close();
    if (outfile.fail()) {
        throw std::runtime_error("Error closing file after writing: " + path);
    }
}

template <BoundaryType BType>
void GridField<BType>::validateValues(double maxAbsValue) {
    for (int i = 0; i < size_; ++i) {
        if (std::abs(values_[i]) > maxAbsValue || std::isnan(values_[i])) {
            throw std::runtime_error("Value validation failed: value out of range or NaN");
        }
    }
}

template <BoundaryType BType>
void GridField<BType>::setRandomInitialCondition(double minValue, double maxValue) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(minValue, maxValue);
    
    for (int i = 0; i < size_; ++i) {
        values_[i] = dist(gen);
    }
}

template <BoundaryType BType>
void GridField<BType>::setLinearInitialCondition() {
    for (int i = 0; i < size_; ++i) {
        // Linear initial condition from -1.0 to 1.0
        values_[i] = -1.0 + 2.0 * i / (size_ - 1.0);
    }
}

template <BoundaryType BType>
void GridField<BType>::setStepInitialCondition() {
    for (int i = 0; i < size_ / 2; ++i) {
        values_[i] = -1.0;
    }
    for (int i = size_ / 2; i < size_; ++i) {
        values_[i] = 1.0;
    }
}

template <BoundaryType BType>
void GridField<BType>::loadInitialConditionWithOffset(const double& offsetValue, const std::string& filePath) {
    std::ifstream file(filePath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open initial condition file: " + filePath);
    }
    
    int tmpSize;
    file.read(reinterpret_cast<char*>(&tmpSize), sizeof(tmpSize));
    if (file.fail()) {
        throw std::runtime_error("Failed to read size from initial condition file");
    }
    
    std::vector<double> tempData(tmpSize);
    file.read(reinterpret_cast<char*>(tempData.data()), tmpSize * sizeof(double));
    if (file.fail()) {
        throw std::runtime_error("Failed to read data from initial condition file");
    }
    
    // Apply read data up to available size
    for (int i = 0; i < std::min(tmpSize, size_); ++i) {
        values_[i] = offsetValue + tempData[i];
    }
    
    // Fill the rest with average value if our grid is larger
    for (int i = tmpSize; i < size_; ++i) {
        values_[i] = offsetValue;
    }
}

template <BoundaryType BType>
void GridField<BType>::loadInitialCondition(const std::string& filePath) {
    std::ifstream file(filePath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open initial condition file: " + filePath);
    }

    int tmpSize;
    file.read(reinterpret_cast<char*>(&tmpSize), sizeof(tmpSize));
    if (file.fail()) {
        throw std::runtime_error("Failed to read size from initial condition file");
    }

    std::vector<double> tempData(tmpSize);
    file.read(reinterpret_cast<char*>(tempData.data()), tmpSize * sizeof(double));
    if (file.fail()) {
        throw std::runtime_error("Failed to read data from initial condition file");
    }

    for (int i = 0; i < size_; ++i) {
        values_[i] = tempData[i];
    }
}

template <BoundaryType BType>
void GridField<BType>::setSinusoidalInitialCondition(double waveNumber, double amplitude, double gridSpacing, double phase) {
    double currentPhase = phase;
    for (int i = 0; i < size_; ++i) {
        values_[i] = amplitude * std::sin(currentPhase);
        currentPhase += waveNumber * gridSpacing;
    }
    
    // Check if the condition is periodic
    double endPhase = phase + waveNumber * gridSpacing * size_;
    if (std::abs(std::sin(phase) - std::sin(endPhase)) > 1e-10) {
        throw std::runtime_error("Initial condition is not periodic");
    }
}

template <BoundaryType BType>
void GridField<BType>::setBookInitialCondition(double domainLength, double gridSpacing) {
    double phi0 = -0.275;
    double amplitude = 0.5;
    
    for (int i = 0; i <= size_ / 2; ++i) {
        values_[i] = phi0 + 2 * amplitude / domainLength * i * gridSpacing - amplitude / 2;
    }
    for (int i = size_ / 2 + 1; i < size_; ++i) {
        values_[i] = phi0 - 2 * amplitude / domainLength * i * gridSpacing + 3 * amplitude / 2;
    }
}

// Явная инстанциация шаблонов для используемых типов граничных условий
template class GridField<BoundaryType::Periodic>;
template class GridField<BoundaryType::Fix>;

// Также нужно явно инстанцировать функцию swap
template void swap(GridField<BoundaryType::Periodic>& a, GridField<BoundaryType::Periodic>& b);
template void swap(GridField<BoundaryType::Fix>& a, GridField<BoundaryType::Fix>& b);

#endif // GRID_FIELD_CPP
