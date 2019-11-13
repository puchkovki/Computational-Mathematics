// Definition of Vector class. That is the class for standard multidimension linear algebra vector
#pragma once

#include <vector>       // std::vector<double>

using namespace std;

// Standard multidimension linear algebra vector
typedef vector<double> Vector;

// Summarizes two vectors with same dimension
Vector operator+(const Vector& left, const Vector& right);

// Substracts two vectors with same dimension
Vector operator-(const Vector& left, const Vector& right);

// Multiplicates two vectors with same dimension by scalar multiplication
double operator*(const Vector& left, const Vector& right);

// Multiplicates vector by number
Vector operator*(const Vector& left, const double right);

// Multiplicates vector by number
Vector operator*(const double left, const Vector& right);

// Divides vector by number
Vector operator/(const Vector& left, const double right);

// Prints vector to output stream os
ostream& operator<<(std::ostream& os, const Vector& to_print);