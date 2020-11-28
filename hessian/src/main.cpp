// Copyright [2019] <Puchkov Kyryll>
#include <iostream>     // Input/output operations
#include <cmath>        // Square root
#include <vector>       // Storing information in std::vector
#include <fstream>      // File managing

// Grid fineness for solution
#define FINENESS 1000

// Table cell width
#define WDTH 20

using namespace std;

// Constant from the task
const double A = 0.4;
// Constant from the task
const double B = 0.2;

// Type for considered functions
typedef double func(double, double);

// Function from the task
double F(double x, double y) {
    double a = 1 - (3 * B) + (3 * A);
    double b = (3 * B) - (6 * A);
    double c = 3 * A;
    double d = sqrt(x * x + y * y) > 1 ? -1 : -sqrt(x * x + y * y);

    vector<double> t(7);
    t[0] = 1.0 / 2.0;
    for (size_t i = 1; i < t.size(); ++i) {
        t[i] = t[i - 1] -
                (a * t[i - 1] * t[i - 1] * t[i - 1] + b * t[i - 1] * t[i - 1] + c * t[i - 1] + d) /
                    (3 * a * t[i - 1] * t[i - 1] + 2 * b * t[i - 1] + c);
    }
    return (1 - t[6]) * (1 - t[6]) * (1 - t[6]) + 3 * (1 - t[6]) * (1 - t[6]) * t[6];
}

// Calculates second mixed derivative
double SecondDerivativeMixed(func f, double x, double y, double h_x, double h_y) {
    return (f(x + h_x, y + h_y) - f(x + h_x, y - h_y) - f(x - h_x, y + h_y) + f(x - h_x, y - h_y)) / (4 * h_x * h_y);
}

// Calculates second derivative by the first variable
double SecondDerivativeX(func f, double x, double y, double h_x) {
    return (f(x + h_x, y) - 2 * f(x, y) + f(x - h_x, y)) / (4 * h_x * h_x);
}

// Calculates second derivative by the second variable
double SecondDerivativeY(func f, double x, double y, double h_y) {
    return (f(x, y + h_y) - 2 * f(x, y) + f(x, y - h_y)) / (4 * h_y * h_y);
}

// Calculates function hessian
double Hessian(func f, double x, double y, double h_x, double h_y) {
    double SDX = SecondDerivativeX(f, x, y, h_x);
    double SDY = SecondDerivativeY(f, x, y, h_y);
    double SDM = SecondDerivativeMixed(f, x, y, h_x, h_y);

    return SDX * SDY - SDM * SDM / 4;
}

int main(void) {
    // Grid cell width
    double step = 2.0 / FINENESS;

    // Working space, storing value of Hessian
    vector<vector<double> > field(FINENESS, vector<double>(FINENESS));

    // Calculating process for each grid node
    for (size_t i = 0; i < FINENESS; ++i) {
        for (size_t j = 0; j < FINENESS; ++j) {
            field[i][j] = Hessian(F, -1.0 + j * step, -1 + i * step, step, step);
        }
    }

    // If Hessian is positive, writes "+", if it's negative, writes "-", else writes "0"
    for (size_t i = 0; i < FINENESS; ++i) {
        for (size_t j = 0; j < FINENESS; ++j) {
            cout << (field[FINENESS - 1 - i][j] > (double) 0 ? "+" : (field[FINENESS - 1 - i][j] < (double) 0 ? "-" : "0"));
        }
        cout << endl;
    }

    /* Values analysis */

    /* TXT file */

    // File output stream
    ofstream values("values.txt");

    // Header
    values.width(WDTH);
    values << "step";
    values.width(WDTH);
    values << "x";
    values.width(WDTH);
    values << "y";
    values.width(WDTH);
    values << "F(x, y)";
    values.width(WDTH);
    values << "F_xx(x, y)";
    values.width(WDTH);
    values << "F_xy(x, y)";
    values.width(WDTH);
    values << "F_yy(x, y)";
    values.width(WDTH);
    values << "Hessian" << endl;

    for (size_t i = 0; i < 8 * WDTH; ++i) {
        values << "-";
    }
    values << endl;

    values.precision(10);
    // Writing
    for (double x = -1; x < -1.0 / 100; x /= 2) {
        for (double y = -1; y < -1.0 / 100; y /= 2) {
            for (size_t i = 10; i <= 100000; i *= 10) {
                step = 2.0 / i;
                values.width(WDTH);
                values << step;
                values.width(WDTH);
                values << x;
                values.width(WDTH);
                values << y;
                values.width(WDTH);
                values << F(x, y);
                values.width(WDTH);
                values << SecondDerivativeX(F, x, y, step);
                values.width(WDTH);
                values << SecondDerivativeMixed(F, x, y, step, step);
                values.width(WDTH);
                values << SecondDerivativeY(F, x, y, step);
                values.width(WDTH);
                values << Hessian(F, x, y, step, step) << endl;
            }
            for (size_t i = 0; i < 8 * WDTH; ++i) {
                values << "-";
            }
            values << endl;
        }
        for (size_t i = 0; i < 8 * WDTH; ++i) {
            values << "-";
        }
        values << endl;
    }

    values.close();

    /* CSV file */

    // File output stream
    values.open("values.csv");

    // Header
    values << "step;";
    values << "x;";
    values << "y;";
    values << "F(x, y);";
    values << "F_xx(x, y);";
    values << "F_xy(x, y);";
    values << "F_yy(x, y);";
    values << "Hessian" << endl;

    values.precision(10);

    // Changing dots to commas
    locale mylocale("");
    values.imbue(mylocale);

    // Writing
    for (double x = -1; x < -1.0 / 100; x /= 2) {
        for (double y = -1; y < -1.0 / 100; y /= 2) {
            for (size_t i = 10; i <= 100000; i *= 10) {
                step = 2.0 / i;
                values << step << ";";
                values << x << ";";
                values << y << ";";
                values << F(x, y) << ";";
                values << SecondDerivativeX(F, x, y, step) << ";";
                values << SecondDerivativeMixed(F, x, y, step, step) << ";";
                values << SecondDerivativeY(F, x, y, step) << ";";
                values << Hessian(F, x, y, step, step) << endl;
            }
        }
    }

    values.close();
}
