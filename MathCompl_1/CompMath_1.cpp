#include <iostream>
#include <cmath>
#include <vector>

#define COLUMN_WIDTH 28

template <typename T>
T Euler_function(T x) {
    T sum = 1;
    for(T I = 1, term = 1; sum != sum + term; sum += term, I++) {
        term *=x / I;
    }
    return sum;
}

template <typename T>
inline T abs(T x) {
    return x > 0 ? x : -x;
}

template <typename T>
T Euler_function_pro(T x) {
    T sum = 1;
    bool inverse = x < 0;
    x = abs(x);
    for(T I = 1, term = 1; sum != sum + term; sum += term, I++) {
        term *=x / I;
    }
    if(inverse) {
        sum = 1 / sum;
    }
    return sum;
}

int main(void) {
    std::vector<double> array = {1, 5, 10, 15, 20, 25, -1, -5, -10, -15, -20, -25};
    
    std::cout.precision(16);
    std::cout.setf(std::ios::left);

    std::cout.width(COLUMN_WIDTH);
    std::cout << "Variable x";

    std::cout.width(COLUMN_WIDTH);
    std::cout << "exp(float)";    

    std::cout.width(COLUMN_WIDTH);
    std::cout << "exp(double)";

    std::cout.width(COLUMN_WIDTH);
    std::cout << "pro_exp(float)";

    std::cout.width(COLUMN_WIDTH);
    std::cout << "pro_exp(double)";

    std::cout.width(COLUMN_WIDTH);
    std::cout << "library exp" << std::endl << std::endl;

    for (auto x: array) {
        std::cout.width(COLUMN_WIDTH);
        std::cout << x;

        //Вычисление значения экспоненты по формуле Тейлора для числа с одинарной точностью (float)
        std::cout.width(COLUMN_WIDTH);
        std::cout  << Euler_function((float)x);

        //Вычисление значения экспоненты по формуле Тейлора для числа с двойной точностью (double)
        std::cout.width(COLUMN_WIDTH);
        std::cout  << Euler_function(x);

        //Вычисление усовершенствованного значения экспоненты по формуле Тейлора для числа с одинарной точностью (float)
        std::cout.width(COLUMN_WIDTH);
        std::cout  << Euler_function_pro((float)x);

        //Вычисление усовершенствованного значения экспоненты по формуле Тейлора для числа с двойной точностью (double)
        std::cout.width(COLUMN_WIDTH);
        std::cout  << Euler_function_pro(x);
        
        //Вычисление табличного значения экспоненты
        std::cout.width(COLUMN_WIDTH);
        std::cout << exp(x) << std::endl << std::endl;
    }
    std::cout << std::endl << std::endl;

    std::cout.width(COLUMN_WIDTH);
    std::cout << "Variable x";
    
    std::cout.width(COLUMN_WIDTH);
    std::cout << "relative error(float)";

    std::cout.width(COLUMN_WIDTH);
    std::cout << "relative error(double)";    

    std::cout.width(COLUMN_WIDTH);
    std::cout << "relative pro_error(float)";
    
    std::cout.width(COLUMN_WIDTH);
    std::cout << "relative pro_error(double)" << std::endl << std::endl;

    for (auto x: array) {        
        std::cout.width(COLUMN_WIDTH);
        std::cout << x;

        std::cout.width(COLUMN_WIDTH);
        std::cout << abs((double)(Euler_function((float)x) - exp(x))) / exp(x) * 100;
        
        std::cout.width(COLUMN_WIDTH);
        std::cout << abs((double)(Euler_function(x) - exp(x))) / exp(x) * 100;

        std::cout.width(COLUMN_WIDTH);
        std::cout << abs((double)(Euler_function_pro((float)x) - exp(x))) / exp(x) * 100;
        
        std::cout.width(COLUMN_WIDTH);
        std::cout << abs((double)(Euler_function_pro(x) - exp(x))) / exp(x) * 100 << std::endl << std::endl;
    }

    return 0;
}
///  