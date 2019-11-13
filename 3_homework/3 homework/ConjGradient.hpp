// Definition of conjugated gradient method method for SLAE
#pragma once

#include "MatrixCSR.hpp"        // MatrixCSR class - compressed sparce row matrices
#include "Vector.hpp"           // Vector class - standard multidimension linear algebra vector
#include "Norms.hpp"            // Norm type

// Makes first k iterations of conjugated gradient method
Vector ConjGradientK(const MatrixCSR& A, const Vector& b, size_t k);

// Solves SLAE with accuracy eps by norm Norm, using conjugated gradient method
Vector ConjGradient(const MatrixCSR& A, const Vector& b, norm Norm, double eps);

// Calculates the dependence of accuracy on iteration number for conjugated gradient method
Vector ConjGradientGraph(const MatrixCSR& A, const Vector& b, norm Norm, double eps);