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
    outfile.write(reinterpret_cast<const char*>(&sizeX_), sizeof(sizeX_));
    outfile.write(reinterpret_cast<const char*>(&sizeY_), sizeof(sizeY_));

    // Write the array data
    outfile.write(reinterpret_cast<const char*>(values_), sizeX_ * sizeY_ * sizeof(double));
    
    outfile.close();
    if (outfile.fail()) {
        throw std::runtime_error("Error closing file after writing: " + path);
    }
}

template <BoundaryType BType>
void GridField<BType>::validateValues(double maxAbsValue) {
    for (int i = 0; i < sizeX_ * sizeY_; ++i) {
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
    
    for (int i = 0; i < sizeX_ * sizeY_; ++i) {
        values_[i] = dist(gen);
    }
}

template <BoundaryType BType>
void GridField<BType>::setLinearInitialCondition() {
    for (int i = 0; i < sizeX_; ++i) {
        for (int j = 0; j < sizeY_; ++j) {
            // Linear gradient from corners
            double normX = static_cast<double>(i) / (sizeX_ - 1);
            double normY = static_cast<double>(j) / (sizeY_ - 1);
            (*this)(i, j) = -1.0 + 2.0 * (normX + normY) / 2.0;
        }
    }
}

template <BoundaryType BType>
void GridField<BType>::setStepInitialCondition() {
    for (int i = 0; i < sizeX_; ++i) {
        for (int j = 0; j < sizeY_; ++j) {
            if (i < sizeX_ / 2 && j < sizeY_ / 2) {
                (*this)(i, j) = -1.0;
            } else {
                (*this)(i, j) = 1.0;
            }
        }
    }
}

template <BoundaryType BType>
void GridField<BType>::loadInitialConditionWithOffset(const double& offsetValue, const std::string& filePath) {
    std::ifstream file(filePath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open initial condition file: " + filePath);
    }
    
    int tmpSizeX, tmpSizeY;
    file.read(reinterpret_cast<char*>(&tmpSizeX), sizeof(tmpSizeX));
    file.read(reinterpret_cast<char*>(&tmpSizeY), sizeof(tmpSizeY));
    if (file.fail()) {
        throw std::runtime_error("Failed to read size from initial condition file");
    }
    
    std::vector<double> tempData(tmpSizeX * tmpSizeY);
    file.read(reinterpret_cast<char*>(tempData.data()), tmpSizeX * tmpSizeY * sizeof(double));
    if (file.fail()) {
        throw std::runtime_error("Failed to read data from initial condition file");
    }
    
    // Apply read data up to available size
    for (int i = 0; i < std::min(tmpSizeX, sizeX_); ++i) {
        for (int j = 0; j < std::min(tmpSizeY, sizeY_); ++j) {
            (*this)(i, j) = offsetValue + tempData[i * tmpSizeY + j];
        }
    }
    
    // Fill the rest with offset value if our grid is larger
    for (int i = 0; i < sizeX_; ++i) {
        for (int j = 0; j < sizeY_; ++j) {
            if (i >= tmpSizeX || j >= tmpSizeY) {
                (*this)(i, j) = offsetValue;
            }
        }
    }
}

template <BoundaryType BType>
void GridField<BType>::setSinusoidalInitialCondition(double waveNumberX, double waveNumberY, double amplitude, 
                                                   double gridSpacingX, double gridSpacingY, double phaseX, double phaseY) {
    for (int i = 0; i < sizeX_; ++i) {
        for (int j = 0; j < sizeY_; ++j) {
            double x = i * gridSpacingX;
            double y = j * gridSpacingY;
            (*this)(i, j) = amplitude * std::sin(waveNumberX * x + phaseX) * std::sin(waveNumberY * y + phaseY);
        }
    }
    
    // Check if the condition is periodic (for simplicity, we only check the corners)
    if (BType == BoundaryType::Periodic) {
        double x0 = 0.0, y0 = 0.0;
        double x1 = sizeX_ * gridSpacingX, y1 = sizeY_ * gridSpacingY;
        
        double val00 = amplitude * std::sin(waveNumberX * x0 + phaseX) * std::sin(waveNumberY * y0 + phaseY);
        double val10 = amplitude * std::sin(waveNumberX * x1 + phaseX) * std::sin(waveNumberY * y0 + phaseY);
        double val01 = amplitude * std::sin(waveNumberX * x0 + phaseX) * std::sin(waveNumberY * y1 + phaseY);
        double val11 = amplitude * std::sin(waveNumberX * x1 + phaseX) * std::sin(waveNumberY * y1 + phaseY);
        
        if (std::abs(val00 - val10) > 1e-10 || std::abs(val00 - val01) > 1e-10 || std::abs(val00 - val11) > 1e-10) {
            throw std::runtime_error("Initial condition is not periodic in both directions");
        }
    }
}

template <BoundaryType BType>
void GridField<BType>::setBookInitialCondition(double domainLengthX, double domainLengthY, double gridSpacingX, double gridSpacingY) {
    double phi0 = -0.275;
    double amplitude = 0.5;
    
    for (int i = 0; i < sizeX_; ++i) {
        for (int j = 0; j < sizeY_; ++j) {
            double x = i * gridSpacingX;
            double y = j * gridSpacingY;
            
            // Extending 1D book example to 2D
            if (x <= domainLengthX / 2 && y <= domainLengthY / 2) {
                (*this)(i, j) = phi0 + 2 * amplitude / domainLengthX * x - amplitude / 2;
            } else {
                (*this)(i, j) = phi0 - 2 * amplitude / domainLengthX * x + 3 * amplitude / 2;
            }
        }
    }
}

template <BoundaryType BType>
void GridField<BType>::setLoadHalfFromFileInitialCondition(const std::string& filePath, double constValue) {
    // Загрузка данных из файла
    std::ifstream file(filePath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open file: " + filePath);
    }
    
    int tmpSizeX, tmpSizeY;
    file.read(reinterpret_cast<char*>(&tmpSizeX), sizeof(tmpSizeX));
    file.read(reinterpret_cast<char*>(&tmpSizeY), sizeof(tmpSizeY));
    if (file.fail()) {
        throw std::runtime_error("Failed to read size from file");
    }

    if (tmpSizeX != sizeX_ || tmpSizeY != sizeY_) {
        throw std::runtime_error("tmpSizeX != halfX || tmpSizeY != sizeY_");
    }

    // Рассчитаем половину сетки
    int halfX = sizeX_ / 2;

    // Читаем данные
    std::vector<double> tempData(tmpSizeX * tmpSizeY);
    file.read(reinterpret_cast<char*>(tempData.data()), tmpSizeX * tmpSizeY * sizeof(double));
    if (file.fail()) {
        throw std::runtime_error("Failed to read data from file");
    }

    // Заполняем левую половину данными из файла
    for (int i = 0; i < halfX; ++i) {
        for (int j = 0; j < sizeY_; ++j) {
            (*this)(i, j) = tempData[i * tmpSizeY + j];
        }
    }

    // Заполняем правую половину константой
    for (int i = halfX; i < sizeX_; ++i) {
        for (int j = 0; j < sizeY_; ++j) {
            (*this)(i, j) = constValue;
        }
    }
}

// Explicit template instantiations
template class GridField<BoundaryType::Periodic>;
template class GridField<BoundaryType::Fix>;

// Also need to explicitly instantiate the swap function
template void swap(GridField<BoundaryType::Periodic>& a, GridField<BoundaryType::Periodic>& b);
template void swap(GridField<BoundaryType::Fix>& a, GridField<BoundaryType::Fix>& b);

#endif // GRID_FIELD_CPP
