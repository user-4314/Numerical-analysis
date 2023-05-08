#include <iostream>
#include <cstdio>
using namespace std;
const int N = 6;
double X[N] = {-2,-1,0,1,2,3},Y[N] = {-7,-1,-1,-1,5,23};//線上的點
double KnownCoefficient[4] = {1,0,-1,-1};//方程式的係數
double DDT[N][N] = {0},FDT[N][N] = {0};

void Interpolation(double x){
    double P = KnownCoefficient[0];
    for(int i=1;i<4;i++)
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

void ForwardDifferenceTable(double x){
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
    printf("bring x = %.2f by Forward Difference Table then y = %f\n",x,y);
}

void BackwardDifferenceTable(double x){
    double s = (x-X[N-1]) / abs( (X[1]-X[0]) );
    double y = FDT[0][N-1];
    for(int i=1;i<N;i++){
        double P = 1;
        for(int j=0;j<i;j++)
            P = P * (s+j) / (j+1);
        y = y + P*FDT[i][N-i-1];
    }
    printf("bring x = %.2f by Backward Difference Table then y = %f\n",x,y);
}

void NewtonForward(double x){
    int pos = N/2;
    double s = (x-X[pos]) / abs( (X[1]-X[0]) ) -1;
    double y = FDT[0][pos];
    for(int i=1;pos>0 || i<N;i++){
        double P = 1;
        if(i%2 ==0)
            pos--;
        else
            s++;
        for(int j=0;j<i;j++)
            P = P * (s-j) / (j+1);
        y = y + P*FDT[i][pos];
    }
    printf("bring x = %.2f by Newton Forward then y = %f\n",x,y);
}

int main()
{
    double x = -3.7;
    Interpolation(x);
    printf("\n");

    for(double i=8;i<10;i+=0.1)
        LagramgeMethod(i);
    printf("\n");

    x = 2.4;
    DivideDifference(x);
    printf("\n");

    for(double i=8;i<10;i+=0.1)
        ForwardDifferenceTable(i);
    printf("\n");

    for(double i=8;i<10;i+=0.1)
        BackwardDifferenceTable(i);
    printf("\n");

    for(double i=8;i<10;i+=0.1)
        NewtonForward(i);
    printf("\n");
    return 0;
}
