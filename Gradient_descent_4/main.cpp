#include <iostream>
#include <vector> 
#include <functional>
#include <cmath>

#define GOLD 1.6180339887498948482
#define EPS 1e-9

using namespace std;

typedef double func1(double);
typedef double func2(double, double);
typedef vector<double> Vector;
typedef Vector vfunc(double, double);

class Matrix {
    private:
        vector<Vector> comps;

    public:

        Matrix(double a00 = 0, double a01 = 0, double a10 = 0, double a11 = 0) {
            comps = vector<Vector>(2, Vector(2, 0));
            comps[0][0] = a00;
            comps[0][1] = a01;
            comps[1][0] = a10;
            comps[1][1] = a11;
        }

        Matrix(Matrix& other) {
            comps = other.Comps();
        }

        Matrix operator=(Matrix& other) {
            comps = other.Comps();
            return other;
        }

        vector<Vector> Comps() {
            return comps;
        }

        Vector operator[](size_t i) {
            return comps[i];
        }

        friend Vector operator*(Matrix left, Vector right);

        ~Matrix() { }
};

Vector operator*(Matrix left, Vector right) {
    Vector res(2, 0);
    res[0] = left[0][0] * right[0] + left[0][1] * right[1];
    res[1] = left[1][0] * right[0] + left[1][1] * right[1];
    return res;
}

Vector operator-(Vector left, Vector right) {
    Vector res(2, 0);
    res[0] = left[0] - right[0];
    res[1] = left[1] - right[1];
    return res;
}

Vector operator*(double left, Vector right) {
    Vector res(2, 0);
    res[0] = left * right[0];
    res[1] = left * right[1];
    return res;
}

ostream& operator<<(ostream& os, Vector to_print) {
    os << "[ ";
    for (auto e : to_print) {
        os << e << " ";
    }
    os << "]";
    return os;
}

double F(double x, double y) {
    return 3 * x * x - 2 * x * sqrt(y) + y - 8 * x + 8;
}

double F_x(function<func2> F, double x, double y, double h) {
    return (F(x + h, y) - F(x - h, y)) / (2 * h);
}

double F_y(function<func2> F, double x, double y, double h) {
    return (F(x, y + h) - F(x, y - h)) / (2 * h);
}

Vector F_Grad(function<func2> F, double x, double y, double h) {
    Vector res(2, 0);
    res[0] = F_x(F, x, y, h);
    res[1] = F_y(F, x, y, h);
    return res;
}

Matrix F_J(function<func2> F, double x, double y, double h) {
    auto F1 = [F, h] (double x, double y) {
        return F_x(F, x, y, h);
    };
    auto F2 = [F, h] (double x, double y) {
        return F_y(F, x, y, h);
    };
    double a00 = F_x(F1, x, y, h);
    double a01 = F_y(F1, x, y, h);
    double a10 = F_x(F2, x, y, h);
    double a11 = F_y(F2, x, y, h);
    return Matrix(a00, a01, a10, a11);
}

double Norm1(Vector x) {
    return abs(x[0]) + abs(x[1]);
}

double Norm2(Vector x) {
    return sqrt(x[0] * x[0] + x[1] * x[1]);
}

double NormC(Vector x) {
    return abs(x[0]) > abs(x[1]) ? abs(x[0]) : abs(x[1]);
}

double min_Gold(function<func1> F, double left, double right, double eps) {
    double a = left, b = right;
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    for (; b - a > eps;) {
        x1 = b - (b - a) / GOLD;
        x2 = a + (b - a) / GOLD;
        y1 = F(x1);
        y2 = F(x2);
        if (y1 >= y2) {
            a = x1;
        }
        else {
            b = x2;
        }
    }
    return y1;
}

double argmin_Gold(function<func1> F, double left, double right, double eps) {
    double a = left, b = right;
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    for (; b - a > eps;) {
        x1 = b - (b - a) / GOLD;
        x2 = a + (b - a) / GOLD;
        y1 = F(x1);
        y2 = F(x2);
        if (y1 >= y2) {
            a = x1;
        }
        else {
            b = x2;
        }
    }
    return x1;
}

Vector min_Great_descent(function<func2> F, double x_0, double y_0, double eps){
    Vector point(2, 0);
    point[0] = x_0;
    point[1] = y_0;
    for (Vector n = F_Grad(F, point[0], point[1], eps * 1e3); Norm2(n) > eps; n = F_Grad(F, point[0], point[1], eps * 1e3)) {
        auto F1 = [F, eps, point] (double tau) {
            Vector grad = F_Grad(F, point[0], point[1], eps);
            return F(point[0] - tau * grad[0], point[1] - tau * grad[1]);
        };
        double tau = argmin_Gold(F1, 0, 1, eps);
        point = point - tau*n;
    }
    Vector res(4, 0);
    res[0] = point[0];
    res[1] = point[1];
    res[2] = F(res[0], res[1]);
    res[3] = Norm2(F_Grad(F, point[0], point[1], eps * 1e3));

    return res;
}

Vector minimal_residual_algorithm(function<func2> F, double x_0, double y_0, double eps){
    Vector point(2, 0);
    point[0] = x_0;
    point[1] = y_0;
    for (Vector n = F_Grad(F, point[0], point[1], eps * 1e3); Norm2(n) > eps; n = F_Grad(F, point[0], point[1], eps * 1e3)) {
        auto F1 = [F, eps, point] (double tau) {
            Vector grad = F_Grad(F, point[0], point[1], eps);
            return Norm2(F_Grad(F, point[0] - tau * grad[0], point[1] - tau * grad[1], eps));
        };
        double tau = argmin_Gold(F1, 0, 1, eps);
        point = point - tau*n;
    }
    Vector res(4, 0);
    res[0] = point[0];
    res[1] = point[1];
    res[2] = F(res[0], res[1]);
    res[3] = Norm2(F_Grad(F, point[0], point[1], eps * 1e3));

    return res;
}

Vector min_Const_tau(function<func2> F, double x_0, double y_0, double tau, double eps) {
    Vector point(2, 0);
    point[0] = x_0;
    point[1] = y_0;
    for (Vector n = F_Grad(F, point[0], point[1], eps * 1e3); Norm2(n) > eps; n = F_Grad(F, point[0], point[1], eps * 1e3)) {
        point = point - tau * n;
    }
    Vector res(4, 0);
    res[0] = point[0];
    res[1] = point[1];
    res[2] = F(res[0], res[1]);
    res[3] = Norm2(F_Grad(F, point[0], point[1], eps * 1e3));

    return res;
}

int main(void) {
    double x_0 = 2.1, y_0 = 4.1;

    Vector res = min_Const_tau(F, x_0, y_0, 0.1, EPS);
    cout << res << endl;

    res = min_Great_descent(F, x_0, y_0, EPS);
    cout << res << endl;

    res = minimal_residual_algorithm(F, x_0, y_0, EPS);
    cout << res << endl;

    return 0;
}