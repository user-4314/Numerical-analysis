#include <iostream>
#include <cstdio>
using namespace std;
const int N = 3;
double X[N] = {-2,0,1},Y[N] = {9,-1,-3};//線上的點
double KnownCoefficient[N] = {1,-3,-1};//方程式的係數

void Interpolation(double x){
    double P = KnownCoefficient[0];
    for(int i=1;i<N;i++)
        P = P * x + KnownCoefficient[i];

    printf("bring %.2f in f(x) by Interpolation y is %f\n",x,P);
}

void LagramgeMethod(double x){
    double sum = 0;
    for(int k=0;k<N;k++){
        double P = 1;
        for(int i=0;i<N;i++){
            if(i!=k)
                P = P * ( (x-X[i])/(X[k]-X[i]) );
            else
                continue;
        }
        sum = sum + P*Y[k];
    }
    printf("bring %.2f in f(x) by Lagramge y is %f\n",x,sum);
}

int main()
{
    double x = -2.7;
    Interpolation(x);

    for(double i=-3;i<3;i+=0.1)
        LagramgeMethod(i);
    return 0;
}
