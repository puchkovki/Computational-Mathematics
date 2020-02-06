#include <iostream>
#include <cmath>

#define EPS 1e-6

using namespace std;

//Определение функции
double Func(double x) {
    return x * x - exp(x);
}

double f_1(double x){
    return x;
}

double f_2(double x){
    return 1;
}

double F_1(double x){
    return f_1(x) * Func(x);
}

double F_2(double x){
    return f_2(x) * Func(x);
}

//Формула Котеса рассчета определенного интеграла для равномерной сетки
double Integral(double a, double b, double h) { //Dot product
    double I = (Func(a) + Func(b)) / 2;

    for(double x_0 = a + h; x_0 < b; x_0 += h) {
        I += Func(x_0);
    }

    return I * h;
}

double Trapezoidal_rule (double h, double a, double b, double delta){
    
    //Находим необходимую мелкость разбиения h для требуемого в условии EPS
    while (delta > EPS / 10) {
        I_0 = Integral(a, b, h);
        I_1 = Integral(a, b, h / 2);

        //Правило Рунге оценки погрешности численных методов
        delta = abs(I_1 - I_0) / 3;

        h /= 2;
    }
}

int main(void) {
    cout << "Integrate f(x) = log(x) / (1 + exp(x)) from 0 to 1:" << endl;

    double a = 1e-15, b = 1, h = 1e-6, delta = 1, I_0 = 0, I_1 = 0;

    

    cout.precision(10);
    cout << I_1 << endl << "delta: " << delta << endl << "WolframAlpha: -0.4387473" << endl;
}