// Definition of Seidel method for SLAE
#pragma once

#include "MatrixCSR.hpp"        // MatrixCSR class - compressed sparce row matrices
#include "Vector.hpp"           // Vector class - standard multidimension linear algebra vector
#include "Norms.hpp"            // Norm type

// Makes first k iterations of Seidel method
Vector SeidelK(const MatrixCSR& A, const Vector& b, size_t k);

// Solves SLAE with accuracy eps by norm Norm, using Seidel method
Vector Seidel(const MatrixCSR& A, const Vector& b, norm Norm, double eps);

// Calculates the dependence of accuracy on iteration number for Seidel method
Vector SeidelGraph(const MatrixCSR& A, const Vector& b, norm Norm, double eps);