#include "grid_field.h"

std::ostream& operator<<(std::ostream& os, const BoundaryType& type) {
    switch (type) {
        case BoundaryType::Periodic:
            os << "Periodic";
            break;
        case BoundaryType::DirichletZero:
            os << "DirichletZero";
            break;
        case BoundaryType::NeumannZero:
            os << "NeumannZero";
            break;
        default:
            throw std::invalid_argument("Unknown BoundaryType value");
    }
    return os;
}

template <BoundaryType BType>
void GridField<BType>::saveToFile(const std::string& path) {
    std::ofstream outfile(path, std::ios::binary);
    if (!outfile) {
        throw std::runtime_error("Failed to open file for writing: " + path);
    }

    outfile.write(reinterpret_cast<const char*>(&size_), sizeof(size_));
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
    
    for (int i = 0; i < std::min(tmpSize, size_); ++i) {
        values_[i] = offsetValue + tempData[i];
    }
    
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

    constexpr double compressionFactor = 1.0;
    double backgroundValue = 0.0;
    const double leftCompressedBoundary = 1.0 - compressionFactor;

    for (int i = 0; i < size_; ++i) {
        double xNew = static_cast<double>(i) / (size_ - 1);  // [0,1]

        if (xNew < leftCompressedBoundary) {
            values_[i] = backgroundValue;
        } else {
            double xi = (xNew - leftCompressedBoundary) / compressionFactor; // [0,1]
            xi = std::max(0.0, std::min(1.0, xi));

            double oldPos = xi * (tmpSize - 1);
            int j0 = static_cast<int>(std::floor(oldPos));
            int j1 = std::min(j0 + 1, tmpSize - 1);
            double alpha = oldPos - j0;

            values_[i] = ((1.0 - alpha) * tempData[j0] + alpha * tempData[j1]) / compressionFactor;
        }
    }
}

template <BoundaryType BType>
void GridField<BType>::setSinusoidalInitialCondition(double waveNumber, double amplitude, double gridSpacing, double phase) {
    double currentPhase = phase;
    for (int i = 0; i < size_; ++i) {
        values_[i] = amplitude * std::sin(currentPhase);
        currentPhase += waveNumber * gridSpacing;
    }
    
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
        values_[i] = phi0 + 2.0 * amplitude / domainLength * i * gridSpacing - amplitude / 2.0;
    }
    for (int i = size_ / 2 + 1; i < size_; ++i) {
        values_[i] = phi0 - 2.0 * amplitude / domainLength * i * gridSpacing + 3.0 * amplitude / 2.0;
    }
}

template class GridField<BoundaryType::Periodic>;
template class GridField<BoundaryType::DirichletZero>;
template class GridField<BoundaryType::NeumannZero>;

template void swap(GridField<BoundaryType::Periodic>& a, GridField<BoundaryType::Periodic>& b);
template void swap(GridField<BoundaryType::DirichletZero>& a, GridField<BoundaryType::DirichletZero>& b);
template void swap(GridField<BoundaryType::NeumannZero>& a, GridField<BoundaryType::NeumannZero>& b);
