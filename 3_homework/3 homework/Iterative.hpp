// Definition of simple single-parameter iterative method method for SLAE
#pragma once

#include "MatrixCSR.hpp"        // MatrixCSR class - compressed sparce row matrices
#include "Vector.hpp"           // Vector class - standard multidimension linear algebra vector
#include "Norms.hpp"            // Norm type

// Makes first k iterations of simple single-parameter iterative method
Vector IterativeK(const MatrixCSR& A, const Vector& b, size_t k);

// Solves SLAE with accuracy eps by norm Norm, using simple single-parameter iterative method
Vector Iterative(const MatrixCSR& A, const Vector& b, norm Norm, double eps);

// Calculates the dependence of accuracy on iteration number for simple single-parameter iterative method
Vector IterativeGraph(const MatrixCSR& A, const Vector& b, norm Norm, double eps);