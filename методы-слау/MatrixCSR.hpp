// Definition of MatrixCSR class. That is the class for compressed sparce row matrices
#pragma once

#include <vector>           // std::vector<size_t>
#include "Vector.hpp"       // Vector class - standard multidimension linear algebra vector

using namespace std;

// Type for array of indexes in matrix
typedef vector<size_t> v_size_t;

// Store for maximal non-diagonal element of matrix
struct Max {
    double value;
    size_t i;
    size_t j;
    double sum;
};

// Store for maximal and minimal Eigenvalue of matrix
struct E_Val {
    double max;
    double min;
};

// Compressed sparce row matrix
class MatrixCSR {
    private:

        /* ---------------------------------------------------------------------------------------*/
        
        /* Fields */

        // Matrix dimension
        size_t dim;

        // Number of non-zero elements
        size_t nnz;

        // Non-zero elements of matrix
        Vector AIJ;

        // Row index
        v_size_t IA;

        // Column index
        v_size_t JA;

        // Diagonal elements
        Vector Diag;

        /* ---------------------------------------------------------------------------------------*/
        
        /* Private auxilliary methods */

        // Standard constructor. Is unavailable
        MatrixCSR(void);

        // Get row with number i, counting from zero
        Vector get_Line(size_t i) const;

        // Get maximal non diagonal element
        Max max_Not_Diag_Abs_El(void) const;

    public:

        /* ---------------------------------------------------------------------------------------*/
        
        /* Constructors */

        // Constructor by path to file with matrix
        MatrixCSR(const string path);

        // Constructor by non-zero elements of matrix, row index and column index
        MatrixCSR(const Vector& AIJ, const v_size_t& IA, const v_size_t& JA);

        // Copy constructor
        MatrixCSR(const MatrixCSR& other);

        /* ---------------------------------------------------------------------------------------*/
        
        /* Assign operator */

        // Assign operator overload
        MatrixCSR operator=(const MatrixCSR& other);

        /* ---------------------------------------------------------------------------------------*/
        
        /* Arythmetic operations */

        // Compares matrices for equality
        friend bool operator==(const MatrixCSR& left, const MatrixCSR& right);

        // Compares matrices for inequality
        friend bool operator!=(const MatrixCSR& left, const MatrixCSR& right);

        // Summarizes two matrices with same dimension
        friend MatrixCSR operator+(const MatrixCSR& left, const MatrixCSR& right);

        // Substracts two matrices with same dimension
        friend MatrixCSR operator-(const MatrixCSR& left, const MatrixCSR& right);

        // Multiplicates two matrices with same dimension
        friend MatrixCSR operator*(const MatrixCSR& left, const MatrixCSR& right);
        
        // Multiplicates matrix and vector with same dimension
        friend Vector operator*(const MatrixCSR& left, const Vector& right);

        // Multiplicates matrix and vector with same dimension with Seidel optimisation
        friend Vector Sei_mult(const MatrixCSR& left, const Vector& right, const Vector& d);

        // Multiplicates matrix by number
        friend MatrixCSR operator*(const MatrixCSR& left, const double right);

        // Multiplicates matrix by number
        friend MatrixCSR operator*(const double left, const MatrixCSR& right);

        // Divides matrix by number
        friend MatrixCSR operator/(const MatrixCSR& left, const double right);

        /* ---------------------------------------------------------------------------------------*/
        
        /* Access */

        // Get matrix dimension
        size_t get_Dim(void) const;

        // Get matrix number of non-zero elements
        size_t get_Nnz(void) const;

        // Read access to matrix non-zero elements
        const Vector& get_AIJ(void) const;

        // Read access to matrix row index
        const v_size_t& get_IA(void) const;

        // Read access to matrix column index
        const v_size_t& get_JA(void) const;

        // Get matrix diagonal elements
        const Vector get_Diag(void) const;

        // Read access to matrix rows
        const Vector operator[](size_t i) const;

        // Sets component with indexes i and j to value
        MatrixCSR Set_IJ(double value, size_t i0, size_t j0);

        // Swaps elements with indexes (i1, j1) and (i2, j2)
        MatrixCSR Swap_Elements(size_t i1, size_t j1, size_t i2, size_t j2);

        /* ---------------------------------------------------------------------------------------*/
        
        /* Matrix special methods */

        // Returns transposed matrix
        MatrixCSR Transpose(void) const;

        // Returns inverse matrix
        MatrixCSR Inverse(void) const;

        // Returns identity matrix with the same dimension
        MatrixCSR Same_Dim_E(void) const;

        // Returns matrix spectral radius
        double Spectral_Radius(void) const;

        // Returns matrix spectral radius with accuracy eps
        double Spectral_Radius(double eps) const;

        // Returns maximal Eigenvalue of matrix
        double Max_Eigen_Value(void) const;

        // Returns minimal Eigenvalue of matrix
        double Min_Eigen_Value(void) const;

        // Returns maximal Eigenvalue of matrix with accuracy eps
        double Max_Eigen_Value(double eps) const;

        // Returns minimal Eigenvalue of matrix with accuracy eps
        double Min_Eigen_Value(double eps) const;

        // Returns maximal and minimal Eigenvalue of matrix in structure with accuracy eps
        E_Val Max_and_Min_Eigen_Value(double eps) const;

        // Returns second norm of matrix
        double M_Norm2(void) const;

        // Returns matrix condition number
        double Cond_Number(void) const;

        // Returns all matrix Eigenvalues
        Vector Jacobi_Eigen_Values() const;

        // Returns all matrix Eigenvalues with accuracy eps
        Vector Jacobi_Eigen_Values(double eps) const;

        // Returns matrix of iteration for Jacobi method
        MatrixCSR B_Jac(void) const;

        /* ---------------------------------------------------------------------------------------*/
        
        /* Printing */

        // Prints matrix to output stream os
        friend ostream& operator<<(std::ostream& os, const MatrixCSR& to_print);

        /* ---------------------------------------------------------------------------------------*/
        
        /* Dectructors */

        // Standard destructor
        ~MatrixCSR();

        /* ---------------------------------------------------------------------------------------*/
};