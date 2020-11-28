// Implementations of common vector norms
#include <cmath>            // sqrt operation

#include "Vector.hpp"       // Vector class - standard multidimension linear algebra vector
#include "Norms.hpp"        // Own header

// First vector norm
double Norm1(const Vector& x) {

    double sum = 0.0;
    for (auto e : x) {
        sum += abs(e);
    }

    return sum;
}

// Second vector norm
double Norm2(const Vector& x) {

    double sum = 0.0;
    for (auto e : x) {
        sum += e * e;
    }

    return sqrt(sum);
}

// Continuous vector norm
double NormC(const Vector& x) {

    double max = 0.0;
    for (auto e : x) {
        if (abs(e) > max) {
            max = abs(e);
        }
    }

    return max;
}