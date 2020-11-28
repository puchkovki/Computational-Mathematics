// Copyright [2019] <Puchkov Kyryll>
#include <iostream>
#include <cmath>

#define EPS 1e-6

using namespace std;

// Определение функции
double Func(double x) {
    return log(x) / (1 + exp(x));
}

// Формула Котеса рассчета определенного интеграла для равномерной сетки
double Integral(double a, double b, double h) {
    double I = (Func(a) + Func(b)) / 2;

    for (double x_0 = a + h; x_0 < b; x_0 += h) {
        I += Func(x_0);
    }

    return I * h;
}

int main(void) {
    cout << "Integrate f(x) = log(x) / (1 + exp(x)) from 0 to 1:" << endl;

    double a = 1e-15, b = 1, h = 1e-6, delta = 1, I_0 = 0, I_1 = 0;

    // Находим необходимую мелкость разбиения h для требуемого в условии EPS
    while (delta > EPS / 10) {
        I_0 = Integral(a, b, h);
        I_1 = Integral(a, b, h / 2);

        // Правило Рунге оценки погрешности численных методов
        delta = abs(I_1 - I_0) / 3;

        h /= 2;
    }

    cout.precision(10);
    cout << I_1 << endl << "delta: " << delta << endl << "WolframAlpha: -0.4387473" << endl;
}
