#include <iostream>
#include <cstdio>
using namespace std;
const int N = 3;
double X[N] = {-2,0,2},Y[N] = {9,-1,-3};//線上的點
double KnownCoefficient[N] = {1,-3,-1};//方程式的係數
double DDT[N][N] = {0},FDT[N][N] = {0};

void Interpolation(double x){
    double P = KnownCoefficient[0];
    for(int i=1;i<N;i++)
        P = P * x + KnownCoefficient[i];

    printf("bring x = %.2f by Interpolation then y = %f\n",x,P);
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
    printf("bring x = %.2f by Lagramge then y = %f\n",x,sum);
}

void DivideDifference(double x){
    int dis = 1;
    for(int i=0;i<N;i++)
        DDT[0][i] = Y[i];
    for(int i=1;i<N;i++){
        for(int j=0;j<N-i;j++)
            DDT[i][j]= (DDT[i-1][j+1]-DDT[i-1][j]) / (X[j+dis]-X[j]);
        dis++;
    }

    double y = DDT[0][0];
    for(int i=1;i<N;i++){
        double P = 1;
        for(int j=0;j<i;j++)
            P = P * (x-X[j]);
        y = y + P*DDT[i][0];
    }
    printf("bring x = %.2f by Divide Difference y = %f\n",x,y);
}


void DifferenceTable(double x){
    for(int i=0;i<N;i++)
        FDT[0][i] = Y[i];
    for(int i=1;i<N;i++){
        for(int j=0;j<N-i;j++)
            FDT[i][j]= FDT[i-1][j+1]-FDT[i-1][j];
    }

    double s = (x-X[0]) / abs( (X[1]-X[0]) );
    double y = FDT[0][0];
    for(int i=1;i<N;i++){
        double P = 1;
        for(int j=0;j<i;j++)
            P = P * (s-j) / (j+1);
        y = y + P*FDT[i][0];
    }
    printf("bring x = %.2f by Front Difference Table then y = %f\n",x,y);
}

int main()
{
    double x = -2.7;
    Interpolation(x);
    printf("\n");

    for(double i=-3;i<3;i+=0.1)
        LagramgeMethod(i);
    printf("\n");

    x = 2.4;
    DivideDifference(x);
    printf("\n");

     x = -0.8;
    DifferenceTable(x);
    printf("\n");
    return 0;
}
