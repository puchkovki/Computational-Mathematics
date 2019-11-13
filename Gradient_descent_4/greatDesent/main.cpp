#include <iostream>
#include "GradientSp.h"

using namespace std;

int main()
{
    double x,y,epsilon;
    std::cout <<"x\n";
    std::cin>>x;
    std::cout <<"y\n";
    std::cin>>y;
    std::cout <<"epsilon\n";
    std::cin>>epsilon;
    GreatDescent(x,y,epsilon);
    return 0;
}

