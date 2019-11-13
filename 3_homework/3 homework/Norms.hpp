// Definition of common vector norms
#pragma once

#include "Vector.hpp"       // Vector class - standard multidimension linear algebra vector

// Norm function type, gets Vector, returns its norm (double)
typedef double norm(const Vector&);

// First vector norm
double Norm1(const Vector& x);

// Second vector norm
double Norm2(const Vector& x);

// Continuous vector norm
double NormInfty(const Vector& x);