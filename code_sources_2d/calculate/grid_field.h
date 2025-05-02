#ifndef GRID_FIELD_H
#define GRID_FIELD_H

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>

enum class BoundaryType {
    Periodic,
    Fix,
};

std::ostream& operator<<(std::ostream& os, const BoundaryType& type);

template <BoundaryType BType = BoundaryType::Periodic>
class GridField {
private:
    int sizeX_ = 0;
    int sizeY_ = 0;
    double* values_ = nullptr;
    double boundaryValue_ = 0.0;

    // Helper to convert 2D index to 1D
    inline int index(int i, int j) const {
        return i * sizeY_ + j;
    }

    // Templated accessor methods for different boundary types
    template <BoundaryType T = BType>
    typename std::enable_if<T == BoundaryType::Periodic, double&>::type
    accessWithBoundary(int i, int j) {
        // For periodic boundary, wrap around the grid
        int wrappedI = (i + sizeX_) % sizeX_;
        int wrappedJ = (j + sizeY_) % sizeY_;
        return values_[index(wrappedI, wrappedJ)];
    }

    template <BoundaryType T = BType>
    typename std::enable_if<T == BoundaryType::Fix, double&>::type
    accessWithBoundary(int i, int j) {
        // For fixed boundary, return boundary value if outside the grid
        if (i < 0 || i >= sizeX_ || j < 0 || j >= sizeY_) {
            return boundaryValue_;
        }
        return values_[index(i, j)];
    }

public:
    GridField() = default;
    
    GridField(int sizeX, int sizeY, double boundaryValue = 0.0):
    sizeX_(sizeX), sizeY_(sizeY), boundaryValue_(boundaryValue)
    {
        values_ = new double[sizeX * sizeY];
    }

    // Delete copy constructor 
    GridField(const GridField&) = delete;

    // Delete assignment
    GridField& operator=(const GridField&) = delete;

    // Add move constructor
    GridField(GridField&& other) noexcept:
    sizeX_(other.sizeX_), sizeY_(other.sizeY_), values_(other.values_),
    boundaryValue_(other.boundaryValue_)
    {
        other.values_ = nullptr;
        other.sizeX_ = 0;
        other.sizeY_ = 0;
    }

    // Add move assignment
    GridField& operator=(GridField&& other) noexcept {
        if (this != &other) {
            delete[] values_;
            values_ = other.values_;
            sizeX_ = other.sizeX_;
            sizeY_ = other.sizeY_;
            boundaryValue_ = other.boundaryValue_;
            other.values_ = nullptr;
            other.sizeX_ = 0;
            other.sizeY_ = 0;
        }
        return *this;
    }

    ~GridField() {
        delete[] values_;
    }

    // Access operator using template specialization
    inline double& operator()(int i, int j) {
        if (i >= 0 && i < sizeX_ && j >= 0 && j < sizeY_) {
            return values_[index(i, j)];
        }
        return accessWithBoundary<BType>(i, j);
    }

    inline double laplacian(int i, int j, double hxSquared, double hySquared) {
        // 2D Laplacian (5-point stencil)
        double dxx = ((*this)(i-1, j) - 2.0 * (*this)(i, j) + (*this)(i+1, j)) / hxSquared;
        double dyy = ((*this)(i, j-1) - 2.0 * (*this)(i, j) + (*this)(i, j+1)) / hySquared;
        return dxx + dyy;
    }

    void swap(GridField& other) {
        std::swap(sizeX_, other.sizeX_);
        std::swap(sizeY_, other.sizeY_);
        std::swap(values_, other.values_);
        std::swap(boundaryValue_, other.boundaryValue_);
    }

    void saveToFile(const std::string& path);
    void validateValues(double maxAbsValue);
    void setRandomInitialCondition(double minValue, double maxValue);
    void setLinearInitialCondition();
    void setStepInitialCondition();
    void loadInitialConditionWithOffset(const double& offsetValue, const std::string& filePath);

    // Methods requiring grid parameters
    void setSinusoidalInitialCondition(double waveNumberX, double waveNumberY, double amplitude, double gridSpacingX, double gridSpacingY, double phaseX = 0.0, double phaseY = 0.0);
    void setBookInitialCondition(double domainLengthX, double domainLengthY, double gridSpacingX, double gridSpacingY);
    void setLoadHalfFromFileInitialCondition(const std::string& filePath, double constValue);
};

// Free function for swap following STL convention
template <BoundaryType BType>
void swap(GridField<BType>& a, GridField<BType>& b) {
    a.swap(b);
}

#endif // GRID_FIELD_H
