// Implementation for MatrixCSR class. That is the class for compressed sparce row matrices
#include <iostream>             // I/O operations
#include <fstream>              // File I/O operations
#include <sstream>              // String stream
#include <cmath>                // sqrt operation

#include "General.hpp"          // General macroses
#include "MatrixCSR.hpp"        // MatrixCSR class - compressed sparce row matrices
#include "Norms.hpp"            // Second norm of vector

// Returns sign of number
double sgn(double x) {
    return x < 0 ? -1 : (x > 0 ? 1 : 0);
}

/* ---------------------------------------------------------------------------------------*/

/* Private auxilliary methods */

// Standard constructor. Is unavailable
MatrixCSR::MatrixCSR(void) { }

// Get row with number i, counting from zero
Vector MatrixCSR::get_Line(size_t i) const {

    Vector result(dim, 0.0);
    size_t left = IA[i];
    size_t right = IA[i + 1];

    for (size_t i = left; i < right; ++i) {
        result[JA[i]] = AIJ[i];
    }

    return result;
}

// Get maximal non diagonal element
Max MatrixCSR::max_Not_Diag_Abs_El(void) const {

    Max res;
    res.value = 0.0;
    res.i = -1;
    res.j = -1;
    res.sum = 0.0;

    for (size_t i = 0; i < dim; ++i) {
        Vector line = (*this)[i];
        for (size_t j = 0; j < dim; ++j) {

            if (j != i && abs(line[j]) > abs(res.value)) {
                res.value = line[j];
                res.i = i;
                res.j = j;
            }

            if (j != i) {
                res.sum += line[j];
            }
        }
    }

    return res;
}

/* ---------------------------------------------------------------------------------------*/

/* Constructors */

// Constructor by path to file with matrix
MatrixCSR::MatrixCSR(const string path) {

    ifstream in(path);
    string line;

    bool ia_read = false, ja_read = false, aij_read = false;

    while(getline(in, line)) {
        stringstream str(line);
        string first_word;
        str >> first_word;

        if (first_word == "n") {
            str >> first_word;
            str >> dim;
            str >> first_word >> first_word;
            str >> nnz;
            continue;
        }

        if (first_word == "VECTOR") {
            str >> first_word;

            if (first_word == "IA") {
                ia_read = true;
                ja_read = false;
                aij_read = false;
                continue;
            }

            if (first_word == "JA") {
                ia_read = false;
                ja_read = true;
                aij_read = false;
                continue;
            }

            if (first_word == "AIJ") {
                ia_read = false;
                ja_read = false;
                aij_read = true;
                continue;
            }
        }

        if (ia_read) {
            IA.push_back(stoi(first_word));
            size_t i;

            while (str >> i) {
                IA.push_back(i);
            }
            continue;
        }

        if (ja_read) {
            JA.push_back(stoi(first_word));
            size_t j;

            while (str >> j) {
                JA.push_back(j);
            }
            continue;
        }

        if (aij_read) {
            double a;
            stringstream sa(first_word);
            sa >> a;
            AIJ.push_back(a);

            while (str >> a) {
                AIJ.push_back(a);
            }
            continue;
        }
    }

    for (size_t i = 0; i < dim; ++i) {
        Diag.push_back((*this)[i][i]);
    }
}

// Constructor by non-zero elements of matrix, row index and column index
MatrixCSR::MatrixCSR(const Vector& AIJ, const v_size_t& IA, const v_size_t& JA) {

    this->dim = IA.size() - 1;
    this->nnz = AIJ.size();
    this->AIJ = AIJ;
    this->IA = IA;
    this->JA = JA;

    for (size_t i = 0; i < dim; ++i) {
        Diag.push_back((*this)[i][i]);
    }
}

// Copy constructor
MatrixCSR::MatrixCSR(const MatrixCSR& other) {

    AIJ = other.AIJ;
    IA = other.IA;
    JA = other.JA;
    dim = other.dim;
    nnz = other.nnz;
    Diag = other.Diag;
}

/* ---------------------------------------------------------------------------------------*/

/* Assign operator */

// Assign operator overload
MatrixCSR MatrixCSR::operator=(const MatrixCSR& other) {

    AIJ = other.AIJ;
    IA = other.IA;
    JA = other.JA;
    dim = other.dim;
    nnz = other.nnz;
    Diag = other.Diag;

    return MatrixCSR(other);
}

/* ---------------------------------------------------------------------------------------*/

/* Arithmetic operations */

// Compares matrices for equality
bool operator==(const MatrixCSR& left, const MatrixCSR& right) {

    if (left.dim != right.dim) {
        return false;
    }

    if (left.nnz != right.nnz) {
        return false;
    }

    for (size_t i = 0; i < left.AIJ.size(); ++i) {
        if (abs(left.AIJ[i] - right.AIJ[i]) > EPS) {
            return false;
        }
    }

    for (size_t i = 0; i < left.IA.size(); ++i) {
        if (left.IA[i] != right.IA[i]) {
            return false;
        }
    }

    for (size_t i = 0; i < left.JA.size(); ++i) {
        if (left.JA[i] != right.JA[i]) {
            return false;
        }
    }

    return true;
}

// Compares matrices for inequality
bool operator!=(const MatrixCSR& left, const MatrixCSR& right) {
    return !(left == right);
}

// Summarizes two matrices with same dimension
MatrixCSR operator+(const MatrixCSR& left, const MatrixCSR& right) {

    if (left.dim != right.dim) {
        perror("MatrixCSR: operator+: dimensions do not match");
        exit(1);
    }

    Vector New_AIJ;
    v_size_t New_IA;
    v_size_t New_JA;
    size_t num = 0;

    New_IA.push_back(num);

    for (size_t i = 0; i < left.dim; ++i) {
        Vector line = left[i] + right[i];

        for (size_t j = 0; j < line.size(); ++j) {
            if (abs(line[j]) > EPS) {
                New_AIJ.push_back(line[j]);
                New_JA.push_back(j);
                ++num;
            }
        }

        New_IA.push_back(num);
    }

    return MatrixCSR(New_AIJ, New_IA, New_JA);
}

// Substracts two matrices with same dimension
MatrixCSR operator-(const MatrixCSR& left, const MatrixCSR& right) {

    if (left.dim != right.dim) {
        perror("MatrixCSR: operator+: dimensions do not match");
        exit(1);
    }

    Vector New_AIJ;
    v_size_t New_IA;
    v_size_t New_JA;
    size_t num = 0;

    New_IA.push_back(num);

    for (size_t i = 0; i < left.dim; ++i) {
        Vector line = left[i] - right[i];

        for (size_t j = 0; j < line.size(); ++j) {
            if (abs(line[j]) > EPS) {
                New_AIJ.push_back(line[j]);
                New_JA.push_back(j);
                ++num;
            }
        }

        New_IA.push_back(num);
    }

    return MatrixCSR(New_AIJ, New_IA, New_JA);
}

// Multiplicates two matrices with same dimension
MatrixCSR operator*(const MatrixCSR& left, const MatrixCSR& right) {

    if (left.dim != right.dim) {
        perror("MatrixCSR: operator*: dimensions do not match");
        exit(1);
    }

    MatrixCSR t_right = right.Transpose();
    Vector New_AIJ;
    v_size_t New_IA;
    v_size_t New_JA;
    size_t num = 0;

    New_IA.push_back(0);

    for (size_t i = 0; i < left.dim; ++i) {
        Vector line = left * t_right[i];

        for (size_t j = 0; j < left.dim; ++j) {
            if (abs(line[j]) > EPS) {
                New_AIJ.push_back(line[j]);
                New_JA.push_back(j);
                ++num;
            }
        }

        New_IA.push_back(num);
    }

    return MatrixCSR(New_AIJ, New_IA, New_JA);
}

// Multiplicates matrix and vector with same dimension
Vector operator*(const MatrixCSR& left, const Vector& right) {

    if (left.dim != right.size()) {
        perror("MatrixCSR: operator*: dimensions do not match");
        exit(1);
    }

    Vector result(right.size(), 0.0);
    
    for (size_t i = 0; i < left.dim; ++i) {
        for (size_t j = left.IA[i]; j < left.IA[i + 1]; ++j) {
            result[i] += left.AIJ[j] * right[left.JA[j]];
        }
    }

    return result;
}

// Multiplicates matrix and vector with same dimension with Seidel optimisation
Vector Sei_mult(const MatrixCSR& left, const Vector& right, const Vector& d) {
    if (left.dim != right.size()) {
        perror("MatrixCSR: operator*: dimensions do not match");
        exit(1);
    }

    Vector result = right;
    
    for (size_t i = 0; i < left.dim; ++i) {
        double sub = result[i];

        for (size_t j = left.IA[i]; j < left.IA[i + 1]; ++j) {
            result[i] += left.AIJ[j] * result[left.JA[j]];
        }

        result[i] += d[i] - sub;
    }

    return result;
}

// Multiplicates matrix by number
MatrixCSR operator*(const MatrixCSR& left, const double right) {
    return MatrixCSR(left.AIJ * right, left.IA, left.JA);
}

// Multiplicates matrix by number
MatrixCSR operator*(const double left, const MatrixCSR& right) {
    return MatrixCSR(left * right.AIJ, right.IA, right.JA);
}

// Divides matrix by number
MatrixCSR operator/(const MatrixCSR& left, const double right) {
    return MatrixCSR(left.AIJ / right, left.IA, left.JA);
}

/* ---------------------------------------------------------------------------------------*/

/* Access */

// Get matrix dimension
size_t MatrixCSR::get_Dim(void) const {
    return dim;
}

// Get matrix number of non-zero elements
size_t MatrixCSR::get_Nnz(void) const {
    return nnz;
}

// Read access to matrix non-zero elements
const Vector& MatrixCSR::get_AIJ(void) const {
    return AIJ;
}

// Read access to matrix row index
const v_size_t& MatrixCSR::get_IA(void) const {
    return IA;
}

// Read access to matrix column index
const v_size_t& MatrixCSR::get_JA(void) const {
    return JA;
}

// Get matrix diagonal elements
const Vector MatrixCSR::get_Diag(void) const {
    return Diag;
}

// Read access to matrix rows
const Vector MatrixCSR::operator[](size_t i) const {
    return this->get_Line(i);
}

// Sets component with indexes i and j to value
MatrixCSR MatrixCSR::Set_IJ(double value, size_t i0, size_t j0) {

    Vector New_AIJ;
    v_size_t New_IA;
    v_size_t New_JA;
    size_t num = 0;

    New_IA.push_back(num);

    for (size_t i = 0; i < dim; ++i) {
        Vector line = (*this)[i];

        for (size_t j = 0; j < dim; ++j) {
            if (i == i0 && j == j0) {
                line[j] = value;
            }

            if (abs(line[j]) > EPS) {
                New_AIJ.push_back(line[j]);
                New_JA.push_back(j);
                ++num;
            }
        }

        New_IA.push_back(num);
    }

    *this = MatrixCSR(New_AIJ, New_IA, New_JA);

    if (i0 == j0) {
        Diag[i0] = value;
    }

    return *this;
}

// Swaps elements with indexes (i1, j1) and (i2, j2)
MatrixCSR MatrixCSR::Swap_Elements(size_t i1, size_t j1, size_t i2, size_t j2) {

    if (i1 > i2) {
        swap(i1, i2);
        swap(j1, j2);
    }

    Vector New_AIJ;
    v_size_t New_IA;
    v_size_t New_JA;
    size_t num = 0;

    New_IA.push_back(num);

    for (size_t i = 0; i < dim; ++i) {
        Vector line = (*this)[i];

        for (size_t j = 0; j < dim; ++j) {
            if (i == i1 && j == j1) {
                line[j] = (*this)[i2][j2];
            }

            if (i == i2 && j == j2) {
                line[j] = (*this)[i1][j1];
            }

            if (abs(line[j]) > EPS) {
                New_AIJ.push_back(line[j]);
                New_JA.push_back(j);
                ++num;
            }
        }
        
        New_IA.push_back(num);
    }

    *this = MatrixCSR(New_AIJ, New_IA, New_JA);

    return MatrixCSR(New_AIJ, New_IA, New_JA);
}

/* ---------------------------------------------------------------------------------------*/

/* Matrix special methods */

// Returns transposed matrix
MatrixCSR MatrixCSR::Transpose(void) const {

    Vector New_AIJ;
    v_size_t New_IA;
    v_size_t New_JA;
    size_t num = 0;

    New_IA.push_back(num);

    for (size_t i = 0; i < dim; ++i) {
        Vector row;

        for (size_t j = 0; j < dim; ++j) {
            row.push_back(this->get_Line(j)[i]);
        }

        for (size_t j = 0; j < dim; ++j) {
            if (row[j] != 0) {
                New_AIJ.push_back(row[j]);
                New_JA.push_back(j);
                ++num;
            }
        }

        New_IA.push_back(num);
    }

    return MatrixCSR(New_AIJ, New_IA, New_JA);
}

// Returns inverse matrix
MatrixCSR MatrixCSR::Inverse(void) const {

    MatrixCSR a = *this;
    MatrixCSR E = this->Same_Dim_E();
    v_size_t where(dim, -1);

    for (size_t col = 0, row = 0; (col < dim) && (row < dim); ++col) {
        size_t sel = row;

        for (size_t i = row; i < dim; ++i) {
            if (abs (a[i][col]) > abs (a[sel][col])) {
                sel = i;
            }
        }

        if (abs (a[sel][col]) < EPS) {
            continue;
        }

        for (size_t i = col; i < dim; ++i) {
            a.Swap_Elements(sel, i, row, i);
            E.Swap_Elements(sel, i, row, i);
        }

        where[col] = row;
        
        for (size_t i = 0; i < dim; ++i) {
            if (i != row) {
                double c = a[i][col] / a[row][col];
                for (size_t j = 0; j < dim; ++j) {
                    double new_a = a[i][j] - a[row][j] * c;
                    a.Set_IJ(new_a, i, j);
                    double new_e = E[i][j] - E[row][j] * c;
                    E.Set_IJ(new_e, i, j);
                }
            }
        }

        row++;
    }

    for (size_t i = 0; i < dim; ++i) {
        double c = a[i][i];

        for (size_t j = 0; j < dim; ++j) {
            double new_a = a[i][j] / c;
            a.Set_IJ(new_a < EPS ? 0 : new_a, i, j);
            double new_e = E[i][j] / c;
            E.Set_IJ(new_e < EPS ? 0 : new_e, i, j);
        }
    }
    return E;
}

// Returns identity matrix with the same dimension
MatrixCSR MatrixCSR::Same_Dim_E(void) const {

    Vector E_AIJ(dim, 1.0);
    v_size_t E_IA(dim + 1, 0.0);
    v_size_t E_JA(dim, 0.0);

    for (size_t i = 0; i < E_JA.size(); ++i) {
        E_JA[i] = i;
        E_IA[i] = i;
    }

    E_IA[dim] = dim;

    return MatrixCSR(E_AIJ, E_IA, E_JA);
}

// Returns matrix spectral radius
double MatrixCSR::Spectral_Radius(void) const {

    Vector r_0(dim);
    srand(time(NULL));

    for (size_t i = 0; i < r_0.size(); ++i) {
        r_0[i] = ((double) rand()) * 1 / RAND_MAX;
    }

    double lambda_0 = 0;
    double lambda_1 = r_0 * (*this * r_0) / (r_0 * r_0);

    while (abs(lambda_0 - lambda_1) > EPS) {
        lambda_0 = lambda_1;
        r_0 = *this * r_0 / Norm2(*this * r_0);
        lambda_1 = r_0 * (*this * r_0) / (r_0 * r_0);
    }

    return lambda_1;
}

// Returns matrix spectral radius with accuracy eps
double MatrixCSR::Spectral_Radius(double eps) const {

    Vector r_0(dim);
    srand(time(NULL));

    for (size_t i = 0; i < r_0.size(); ++i) {
        r_0[i] = ((double) rand()) * 10 / RAND_MAX;
    }

    double lambda_0 = 0;
    double lambda_1 = r_0 * (*this * r_0) / (r_0 * r_0);

    while (abs(lambda_0 - lambda_1) > eps) {
        lambda_0 = lambda_1;
        r_0 = *this * r_0 / Norm2(*this * r_0);
        lambda_1 = r_0 * (*this * r_0) / (r_0 * r_0);
    }

    return lambda_1;
}

// Returns maximal Eigenvalue of matrix
double MatrixCSR::Max_Eigen_Value(void) const {

    Vector e_val = this->Jacobi_Eigen_Values();
    double max = -INF;

    for (auto e : e_val) {
        if (e > max) {
            max = e;
        }
    }

    return max;
}

// Returns minimal Eigenvalue of matrix with accuracy eps
double MatrixCSR::Min_Eigen_Value(void) const {

    Vector e_val = this->Jacobi_Eigen_Values();
    double min = INF;

    for (auto e : e_val) {
        if (e < min) {
            min = e;
        }
    }

    return min;
}

// Returns maximal Eigenvalue of matrix with accuracy eps
double MatrixCSR::Max_Eigen_Value(double eps) const {

    Vector e_val = this->Jacobi_Eigen_Values(eps);
    double max = -INF;

    for (auto e : e_val) {
        if (e > max) {
            max = e;
        }
    }

    return max;
}

// Returns minimal Eigenvalue of matrix
double MatrixCSR::Min_Eigen_Value(double eps) const {

    Vector e_val = this->Jacobi_Eigen_Values(eps);
    double min = INF;

    for (auto e : e_val) {
        if (e < min) {
            min = e;
        }
    }

    return min;
}

// Returns maximal and minimal Eigenvalue of matrix in structure  with accuracy eps
E_Val MatrixCSR::Max_and_Min_Eigen_Value(double eps) const {

    E_Val res;
    res.max = this->Spectral_Radius(eps);
    double temp = (res.max * this->Same_Dim_E() - *this).Spectral_Radius(eps);
    res.min = res.max - temp;

    return res;
}

// Returns second norm of matrix
double MatrixCSR::M_Norm2(void) const {

    MatrixCSR a = *this;
    double e_val;

    if (a == a.Transpose()) {
        e_val = a.Max_Eigen_Value();
        return e_val;
    }

    e_val = (a.Transpose() * a).Max_Eigen_Value();

    return sqrt(e_val);
}

// Returns matrix condition number
double MatrixCSR::Cond_Number(void) const {
    return this->M_Norm2() * this->Inverse().M_Norm2();
}

// Returns all matrix Eigenvalues
Vector MatrixCSR::Jacobi_Eigen_Values(void) const {

    MatrixCSR A_0 = *this;

    if (A_0 != A_0.Transpose()) {
        perror("MatrixCSR: Jacobi_Eigen_Values: Matrix is not symmetric");
    }

    Max aij = A_0.max_Not_Diag_Abs_El();

    while (abs(aij.value) > EPS * 100000) {
        //cout << abs(aij.value) << " " << aij.sum << endl;
        double den = A_0[aij.i][aij.i] - A_0[aij.j][aij.j];
        double tan_2_phi = den < EPS ? INF : 2 * aij.value / den;
        double sin_phi = tan_2_phi < INF ? sgn(tan_2_phi) * sqrt(0.5 * (1 - 1 / sqrt(1 + tan_2_phi * tan_2_phi))) : 1 / sqrt(2);
        double cos_phi = tan_2_phi < INF ? sqrt(0.5 * (1 + 1 / sqrt(1 + tan_2_phi * tan_2_phi))) : 1 / sqrt(2);

        Vector H_AIJ;
        v_size_t H_IA;
        v_size_t H_JA;
        size_t num = 0;

        H_IA.push_back(num);

        for (size_t i = 0; i < dim; ++i) {
            if (i == aij.i) {
                H_AIJ.push_back(cos_phi);
                H_JA.push_back(aij.i);
                H_AIJ.push_back(-sin_phi);
                H_JA.push_back(aij.j);
                num += 2;
                H_IA.push_back(num);
                continue;
            }

            if (i == aij.j) {
                H_AIJ.push_back(sin_phi);
                H_JA.push_back(aij.i);
                H_AIJ.push_back(cos_phi);
                H_JA.push_back(aij.j);
                num += 2;
                H_IA.push_back(num);
                continue;
            }

            H_AIJ.push_back(1.0);
            H_JA.push_back(i);
            ++num;
            H_IA.push_back(num);
        }

        MatrixCSR H(H_AIJ, H_IA, H_JA);
        
        A_0 = (H.Transpose() * A_0) * H;
        aij = A_0.max_Not_Diag_Abs_El();
    }

    Vector res;
    for (size_t i = 0; i < dim; ++i) {
        res.push_back(A_0.get_Diag()[i]);
    }

    return res;
}

// Returns all matrix Eigenvalues with accuracy eps
Vector MatrixCSR::Jacobi_Eigen_Values(double eps) const {

    MatrixCSR A_0 = *this;

    if (A_0 != A_0.Transpose()) {
        perror("MatrixCSR: Jacobi_Eigen_Values: Matrix is not symmetric");
    }

    Max aij = A_0.max_Not_Diag_Abs_El();

    while (abs(aij.value) > eps) {
        double den = A_0[aij.i][aij.i] - A_0[aij.j][aij.j];
        double tan_2_phi = den < EPS ? INF : 2 * aij.value / den;
        double sin_phi = tan_2_phi < INF ? sgn(tan_2_phi) * sqrt(0.5 * (1 - 1 / sqrt(1 + tan_2_phi * tan_2_phi))) : 1 / sqrt(2);
        double cos_phi = tan_2_phi < INF ? sqrt(0.5 * (1 + 1 / sqrt(1 + tan_2_phi * tan_2_phi))) : 1 / sqrt(2);

        Vector H_AIJ;
        v_size_t H_IA;
        v_size_t H_JA;
        size_t num = 0;

        H_IA.push_back(num);

        for (size_t i = 0; i < dim; ++i) {
            if (i == aij.i) {
                H_AIJ.push_back(cos_phi);
                H_JA.push_back(aij.i);
                H_AIJ.push_back(-sin_phi);
                H_JA.push_back(aij.j);
                num += 2;
                H_IA.push_back(num);
                continue;
            }

            if (i == aij.j) {
                H_AIJ.push_back(sin_phi);
                H_JA.push_back(aij.i);
                H_AIJ.push_back(cos_phi);
                H_JA.push_back(aij.j);
                num += 2;
                H_IA.push_back(num);
                continue;
            }

            H_AIJ.push_back(1.0);
            H_JA.push_back(i);
            ++num;
            H_IA.push_back(num);
        }

        MatrixCSR H(H_AIJ, H_IA, H_JA);

        A_0 = (H.Transpose() * A_0) * H;
        aij = A_0.max_Not_Diag_Abs_El();
    }

    Vector res;
    for (size_t i = 0; i < dim; ++i) {
        res.push_back(A_0.get_Diag()[i]);
    }

    return res;
}

// Returns matrix of iteration for Jacobi method
MatrixCSR MatrixCSR::B_Jac(void) const {

    MatrixCSR B = *this;
    
    for (size_t i = 0; i < B.dim; ++i) {
        for (size_t j = B.IA[i]; j < B.IA[i + 1]; ++j) {
            if (B.JA[j] == i) {
                B.AIJ[j] = 0;
            }
            else {
                B.AIJ[j] /= -Diag[i];
            }
        }
    }

    B.Diag = Vector(dim, 0.0);

    return B;
}

/* ---------------------------------------------------------------------------------------*/

/* Printing */

// Prints matrix to output stream os
ostream& operator<<(ostream& os, const MatrixCSR& to_print) {

    for (size_t i = 0; i < to_print.dim; ++i) {
        os << to_print[i] << endl;
    }
    
    return os;
}

/* ---------------------------------------------------------------------------------------*/

/* Destructors */

// Standard destructor
MatrixCSR::~MatrixCSR() { }
