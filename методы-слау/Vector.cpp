// Implementation of Vector class. That is the class for standard multidimension linear algebra vector
#include <iostream>         // I/O operations

#include "Vector.hpp"       // Own header

// Summarizes two vectors with same dimension
Vector operator+(const Vector& left, const Vector& right) {

    if (left.size() != right.size()) {
        perror("std::vector<double>: operator+: dimensions do not match");
        exit(1);
    }

    size_t size = left.size();
    Vector result(size, 0.0);

    for (size_t i = 0; i < size; ++i) {
        result[i] = left[i] + right[i];
    }

    return result;
}

// Substracts two vectors with same dimension
Vector operator-(const Vector& left, const Vector& right) {

    if (left.size() != right.size()) {
        perror("std::vector<double>: operator-: dimensions do not match");
        exit(1);
    }

    size_t size = left.size();
    Vector result(size, 0.0);

    for (size_t i = 0; i < size; ++i) {
        result[i] = left[i] - right[i];
    }

    return result;
}

// Multiplicates two vectors with same dimension by scalar multiplication
double operator*(const Vector& left, const Vector& right) {

    if (left.size() != right.size()) {
        perror("std::vector<double>: operator*: dimensions do not match");
        exit(1);
    }

    double sum = 0.0;
    size_t size = left.size();

    for (size_t i = 0; i < size; ++i) {
        sum += left[i] * right[i];
    }

    return sum;
}

// Multiplicates vector by number
Vector operator*(const Vector& left, const double right) {

    size_t size = left.size();
    Vector result(size, 0.0);

    for (size_t i = 0; i < size; ++i) {
        result[i] = left[i] * right;
    }

    return result;
}

// Multiplicates vector by number
Vector operator*(const double left, const Vector& right) {

    size_t size = right.size();
    Vector result(size, 0.0);

    for (size_t i = 0; i < size; ++i) {
        result[i] = left * right[i];
    }

    return result;
}

// Divides vector by number
Vector operator/(const Vector& left, const double right) {

    size_t size = left.size();
    Vector result(size, 0.0);

    for (size_t i = 0; i < size; ++i) {
        result[i] = left[i] / right;
    }

    return result;
}

// Prints vector to output stream os
ostream& operator<<(ostream& os, const Vector& to_print) {

    for (auto e : to_print) {
        os.width(25);
        os << e << " ";
    }
    
    return os;
}