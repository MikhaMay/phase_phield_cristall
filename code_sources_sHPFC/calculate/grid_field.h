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
    int size_ = 0;
    double* values_ = nullptr;
    double leftBoundary_ = 0.0;
    double rightBoundary_ = 0.0;

    // Templated accessor methods for different boundary types
    template <BoundaryType T = BType>
    typename std::enable_if<T == BoundaryType::Periodic, double&>::type
    accessWithBoundary(int i) {
        if (i == -1) return values_[size_ - 1];
        if (i == size_) return values_[0];
        return values_[i];
    }

    template <BoundaryType T = BType>
    typename std::enable_if<T == BoundaryType::Fix, double&>::type
    accessWithBoundary(int i) {
        if (i == -1) return leftBoundary_;
        if (i == size_) return rightBoundary_;
        return values_[i];
    }

    // Templated accessor methods for different boundary types (const version)
    template <BoundaryType T = BType>
    typename std::enable_if<T == BoundaryType::Periodic, const double&>::type
    accessWithBoundary(int i) const {
        if (i == -1) return values_[size_ - 1];
        if (i == size_) return values_[0];
        return values_[i];
    }

    template <BoundaryType T = BType>
    typename std::enable_if<T == BoundaryType::Fix, const double&>::type
    accessWithBoundary(int i) const {
        if (i == -1) return leftBoundary_;
        if (i == size_) return rightBoundary_;
        return values_[i];
    }

public:
    GridField() = default;
    
    GridField(int size, double leftBoundary, double rightBoundary):
    size_(size), leftBoundary_(leftBoundary), rightBoundary_(rightBoundary)
    {
        values_ = new double[size];
    }

    // Delete copy constructor 
    GridField(const GridField&) = delete;

    // Delete assignment
    GridField& operator=(const GridField&) = delete;

    // Add move constructor
    GridField(GridField&& other) noexcept:
    size_(other.size_), values_(other.values_),
    leftBoundary_(other.leftBoundary_), rightBoundary_(other.rightBoundary_)
    {
        other.values_ = nullptr;
        other.size_ = 0;
    }

    // Add move assignment
    GridField& operator=(GridField&& other) noexcept {
        if (this != &other) {
            delete[] values_;
            values_ = other.values_;
            size_ = other.size_;
            leftBoundary_ = other.leftBoundary_;
            rightBoundary_ = other.rightBoundary_;
            other.values_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }

    ~GridField() {
        delete[] values_;
    }

    // Accessor methods
    const int& size() const { return size_; }

    // Bracket operator using template specialization
    inline double& operator[] (int i) {
        if (i >= 0 && i < size_) {
            return values_[i];
        }
        return accessWithBoundary<BType>(i);
    }

    inline const double& operator[] (int i) const {
        if (i >= 0 && i < size_) {
            return values_[i];
        }
        return accessWithBoundary<BType>(i);
    }

    inline double laplacian(int i, double hSquared) const {
        return ((*this)[i-1] - 2 * (*this)[i] + (*this)[i+1]) / hSquared;
    }

    inline double nabla(int i, double h, bool forward = true) const {
        if (forward) {
            return ((*this)[i+1] - (*this)[i]) / h;
        } else {
            return ((*this)[i] - (*this)[i-1]) / h;
        }
    }

    void swap(GridField& other) {
        std::swap(size_, other.size_);
        std::swap(values_, other.values_);
        std::swap(leftBoundary_, other.leftBoundary_);
        std::swap(rightBoundary_, other.rightBoundary_);
    }

    void saveToFile(const std::string& path);
    void validateValues(double maxAbsValue);
    void setRandomInitialCondition(double minValue, double maxValue);
    void setLinearInitialCondition();
    void setStepInitialCondition();
    void loadInitialConditionWithOffset(const double& offsetValue, const std::string& filePath);
    void loadInitialCondition(const std::string& filePath);

    // Методы, требующие параметры сетки, передадим их как аргументы
    void setSinusoidalInitialCondition(double waveNumber, double amplitude, double gridSpacing, double phase = 0.0);
    void setBookInitialCondition(double domainLength, double gridSpacing);
};

// Free function for swap following STL convention
template <BoundaryType BType>
void swap(GridField<BType>& a, GridField<BType>& b) {
    a.swap(b);
}

#endif // GRID_FIELD_H
