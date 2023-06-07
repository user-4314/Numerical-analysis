#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
using namespace std;

const int N = 5;
double X[N] = {-2,0,1,2,3} , Y[N] = {27,7,6,11,22};
int sizeofarray;
double x[N][N+1]={0} , a[N][N]={0} , O[N] = {0};//x[0][n+1]~x[n][n+1] represent y ,a[k][i] represent Pk of coefficient Ai ,O[n] represent the error

void PartialPivoting(int pos){
    //change the biggest to the first
    float ma = x[pos][pos];
    int mapos = pos;
    for(int i = pos+1;i < sizeofarray+1;i++){
        if(abs(x[i][0]) > ma){
            ma = abs(x[i][0]);
            mapos = i;
        }
    }
    for(int j = 0;j < N+1;j++){
        float counter = x[mapos][j];
        x[mapos][j] = x[pos][j];
        x[pos][j] = counter;
    }
}

void SolveLinearEquation(int sizeofarray){
    //after Gaussian elimination
    for(int k = 0;k < sizeofarray;k++){
        PartialPivoting(k);
        for(int i = k+1;i < sizeofarray+1;i++){
            float d = x[i][k]/x[k][k];
            for(int j = k;j <= sizeofarray+1;j++){
                x[i][j] = x[i][j] - d*x[k][j];
            }
        }
    }
    //debug print
    /*for(int i = 0;i < sizeofarray+1;i++){
        for(int j = 0;j < sizeofarray+2;j++)
            printf("%7.2f ",x[i][j]);
        printf("\n");
    }
    printf("\n");*/
    //the solution
    for(int k = sizeofarray;k>=0;k--){
        float sum = 0;
        for(int j = k+1;j < sizeofarray+1;j++){
            sum = sum + x[k][j]*a[sizeofarray][j];
        }
        a[sizeofarray][k] = (x[k][sizeofarray+1] - sum) / x[k][k];
    }
    //debug print
    for(int i = 0;i < sizeofarray+1;i++)
        printf("a%d is %f\n",i,a[sizeofarray][i]);
}

void MadeLeastSquae(int n){
    memset(x,0,sizeof(x));

    if(n <= N-1){
        //made the x array
        for(int i = 0;i <= n;i++){
            for(int j = 0;j <= n;j++){
                //sum the total x
                for(int k = 0;k < N;k++)
                    x[i][j] = x[i][j] + pow(X[k],i+j);
            }
        }
        //made the y array
        for(int i = 0;i <= n;i++){
            for(int k = 0;k < N;k++)
                x[i][n+1] = x[i][n+1] + Y[k] * pow(X[k],i);
        }
        //debug print
        /*for(int i = 0;i <= n;i++){
            for(int j = 0;j <= n+1;j++){
                printf("%7.2f ",x[i][j]);
            }
            printf("\n");
        }*/
        printf("\n");
        //solve the ans of coefficient
        sizeofarray = n;
        SolveLinearEquation(sizeofarray);
    }
    else
        printf("\npoint is not enough to solve at this degree!!\n");
}

double Interpolation(int n,double x){
    //solve¥XPn(x)
    double P = a[n][n];
    for(int i=n-1;i>=0;i--)
        P = P * x + a[n][i];

   return P;
}

void calculateO(int n){
    //big O formula sum(yi-P(xi))^2/n-m
    for(int i=0;i<N;i++)
        O[n] = O[n] + (Y[i] - Interpolation(n,X[i]))*(Y[i] - Interpolation(n,X[i]));

    O[n] = O[n]/ (N-n);
    printf("Error big O is %f\n",O[n]);
    printf("\n");
}


int main()
{
    for(int i = 1;i <= N;i++){
        printf("the P%d coefficient is",i);
        MadeLeastSquae(i);
        calculateO(i);

        //if the next Pi can not provide more precise then stop
        if(O[i-1]/O[i]<2 && i>1){
            printf("P%d is the best choice.\n",i-1);
            break;
        }
        //illustrate the polynomial
        /*for(double x = -2;x<=3.1;x+=0.1)
            printf("P%d(%5.2f) = %f\n",i,x,Interpolation(i,x));
        printf("\n");*/
    }

    int n;
    printf("enter the exponential you want:");
    while(scanf("%d",&n)!=EOF){
        printf("the polynomial P%d is like:",n);
        for(int i = 0;i < n+1;i++){
            if(i == 0)
                printf("%f ",a[n][i]);
            else{
                if(a[n][i]<0)
                    printf(" %fx^%d",a[n][i],i);
                else
                    printf(" + %fx^%d",a[n][i],i);
            }
        }

        printf("\n");
        for(double x = -2;x<=3.1;x+=0.1)
            printf("P%d(%5.2f) = %f\n",n,x,Interpolation(n,x));
        printf("\nenter the exponential you want:");

    }
    return 0;
}
