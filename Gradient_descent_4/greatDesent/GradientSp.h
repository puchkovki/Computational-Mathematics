#ifndef GRADIENTSP_H
#define GRADIENTSP_H

#define NMAX 1000

#include<cmath>
#include<cstdio>
#include<iostream>

//Собственно здесь записывается наша функция
double f(double x,double y)
{
    return x*x + 2*x*y + 3*y*y - 2*x -3*y;
}

//Это первая производная по dx
double f_dx(double x,double y)
{
    return 2*x + 2*y -2;
}

//Это первая производная по dy
double f_dy(double x,double y)
{
    return 2*x + 6*y -3;
}

//Это функция g d методе наискорейшего (градиентного) спуска
double g(double x, double y, double alpha)
{
    return f(x - alpha*f_dx(x,y), y - alpha*f_dy(x,y));
}

//двумерная норма
double norma(double x, double y)
{
    return sqrt(x*x+y*y);
}

//Метод половинного деления для нахождения минимума в градиентном спуске
double Dihotomia(double a0, double b0, double epsilon, double x, double y)
{
    //Номер шага
    int k;    
    //Отклонени от середины отрезка влево, вправо
    double lk, mk;
    //Величина на которую мы отклонимся от середины отрезка
    double delta=0.5*epsilon;
    //Точка минимума
    double x_;
    //Отрезок локализации минимума
    double ak=a0, bk=b0;    
    k=1;
    //Пока длина отрезка больше заданной точности
    do
    {
        //Берем середину (ну почти середину - +\- некоторое дельта в частности у нас delta=0.5*epsilon)
        lk=(ak+bk-delta)/2;
        mk=(ak+bk+delta)/2;

        k++;
        //Проверяем в какую часть попадает точка минимума слева от разбиения или справа и выбираем соответствующую точку
        if(g(x,y,lk)<=g(x,y,mk))
        {
            //Теперь правая граница отрезка локализации равна mk
            bk=mk;
        }
        else
        {
            //Теперь левая граница отрезка локализации равна mk
            ak=lk;
        }
    } while ((bk-ak)>=epsilon);

    x_=(ak+bk)/2; //minimum point

    return x_;
}


// метод наискорейшего спуска
double GreatDescent(int bx, int by,double epsilon)
{
    double x[NMAX];
    double y[NMAX];
    double alpha[NMAX];
    int k;

    //Начальное приближение u[0]
    x[0]=bx;
    y[0]=by;

    std::cout <<"Результаты:"
            <<std::endl<< "x("<<0<<"): ("<<x[0]<<", "<<y[0]<<")"<<std::endl;

    for (k=0; ; k++)
    {
        //Находим alpha_k как минимум функции g на отрезке -10000,100000
        alpha[k]=Dihotomia(-10000,100000,epsilon,x[k],y[k]);
        //Вычисляем u[k]
        x[k+1]=x[k]-alpha[k]*f_dx(x[k], y[k]);
        y[k+1]=y[k]-alpha[k]*f_dy(x[k], y[k]);

        std::cout<<"x("<<k+1<<"): "<<"("<<x[k+1]<<", "<<y[k+1]<<")"<<std::endl
                   <<"f("<<x[k+1]<<", "<<y[k+1]<<") = "<<f(x[k+1], y[k+1]) <<std::endl;
        if (k>1)
        {
            //Проверяем условие остановки
            if (norma(x[k+1]-x[k],y[k+1]-y[k])<epsilon)
            {
                break;
            }
        }
    }

    std::cout <<"Точка минимума (epsilon="<<epsilon<<")"<<std::endl
             <<"f("<<x[k+1]<<", "<<y[k+1]<<") = "<<f(x[k+1], y[k+1]) <<std::endl;
    return f(x[k+1],y[k+1]);
}
#endif // GRADIENTSP_H
