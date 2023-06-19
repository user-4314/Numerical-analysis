#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
using namespace std;

const int N = 43;
double X[N] = {0.00 , 0.14 , 0.45 , 0.73 , 0.88 , 1.05 , 1.33 , 1.64 , 1.92 , 2.06 ,
                2.37 , 2.54 , 2.74 , 2.93 , 3.24 , 3.53 , 3.79 , 4.11 , 4.39 , 4.54 ,
                4.69 , 5.02 , 5.32 , 5.56 , 5.83 , 5.91 , 5.99 , 6.07 , 6.32 , 6.54 ,
                6.73 , 6.96 , 7.22 , 7.31 , 7.46 , 7.68 , 7.82 , 8.12 , 8.22 , 8.43 ,
                8.67 , 8.96 , 9.19 } ,
Y[N] = {1.0000, 1.1728 , 1.6547 , 2.0248 , 2.1168 , 2.0994 , 1.7931 , 1.1986, 0.7258 , 0.5851 ,
        0.6046, 0.7890 , 1.0916 , 1.3965 , 1.7492 , 1.8157 , 1.6898 , 1.4771, 1.3591 , 1.3375 ,
        1.3333, 1.3104 , 1.1899 , 1.0226 , 0.8710 , 0.8488 , 0.8441 , 0.8598, 1.0359 , 1.3468 ,
        1.6390, 1.9657 , 2.1259 , 2.1096 , 1.9918 , 1.6682 , 1.4028 , 0.8473, 0.7065 , 0.5483 ,
        0.6205, 0.9991 , 1.3630 } ;
int sizeofarray;
long double x[N][N+1]={0} , a[N][N]={0} , O[N] = {0};
//x[0][n+1]~x[n][n+1] represent y ,a[k][i] represent Pk of coefficient Ai ,O[n] represent the error

void PartialPivoting(int pos){
    //change the biggest to the first
    double ma = x[pos][pos];
    int mapos = pos;
    for(int i = pos+1;i < sizeofarray+1;i++){
        if(abs(x[i][0]) > ma){
            ma = abs(x[i][0]);
            mapos = i;
        }
    }
    for(int j = 0;j < N+1;j++){
        double counter = x[mapos][j];
        x[mapos][j] = x[pos][j];
        x[pos][j] = counter;
    }
}

void SolveLinearEquation(int sizeofarray){
    //after Gaussian elimination
    for(int k = 0;k < sizeofarray;k++){
        PartialPivoting(k);
        for(int i = k+1;i < sizeofarray+1;i++){
            double d = x[i][k]/x[k][k];
            for(int j = k;j <= sizeofarray+1;j++){
                x[i][j] = x[i][j] - d*x[k][j];
            }
        }
    }
    //debug print
    /*for(int i = 0;i < sizeofarray+1;i++){
        for(int j = 0;j < sizeofarray+2;j++)
            printf("%7.2Lf ",x[i][j]);
        printf("\n");
    }
    printf("\n");*/
    //the solution
    for(int k = sizeofarray;k>=0;k--){
        double sum = 0;
        for(int j = k+1;j < sizeofarray+1;j++){
            sum = sum + x[k][j]*a[sizeofarray][j];
        }
        a[sizeofarray][k] = (x[k][sizeofarray+1] - sum) / x[k][k];
    }
    //debug print
    for(int i = 0;i < sizeofarray+1;i++)
        printf("a%d is %Lf\n",i,a[sizeofarray][i]);
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
                printf("%Lf ",x[i][j]);
            }
            printf("\n");
        }
        printf("\n");*/
        //solve the ans of coefficient
        sizeofarray = n;
        SolveLinearEquation(sizeofarray);
    }
    else
        printf("\npoint is not enough to solve at this degree!!\n");
}

double Interpolation(int n,double x){
    //solve XPn(x)
    double P = a[n][n];
    for(int i=n-1;i>=0;i--)
        P = P * x + a[n][i];

   return P;
}

void calculateO(int n){
    //big O formula sum(yi-P(xi))^2/n-m
    for(int i=0;i<N;i++){
        O[n] = O[n] + (Y[i] - Interpolation(n,X[i]))*(Y[i] - Interpolation(n,X[i]));
        //printf("error sum to x%2d is %Lf\n",i+1,O[n]);
    }

    //O[n] = O[n]/ (N-n);
    printf("Error big O is %Lf\n",O[n]);
    printf("\n");
}


int main()
{
    for(int i = 1;i <= N;i++){
        printf("the P%d coefficient is\n",i);
        MadeLeastSquae(i);
        calculateO(i);

        //if the next Pi can not provide more precise then stop
        /*if(O[i-1]/O[i]<2 && i>1){
            printf("P%d is the best choice.\n",i-1);
            break;
        }*/
    }

    int n;
    printf("enter the exponential you want:");
    //choose the exponential
    while(scanf("%d",&n)!=EOF){
        printf("the polynomial P%d is like:",n);
        //show the entire polynomial
        for(int i = 0;i < n+1;i++){
            if(i == 0)
                printf("%Lf ",a[n][i]);
            else{
                if(a[n][i]<0)
                    printf(" %Lfx^%d",a[n][i],i);
                else
                    printf(" + %Lfx^%d",a[n][i],i);
            }
        }

        printf("\n");
        //illustrate the polynomial
        for(int x = 0;x<N;x++)
            printf("P%d(%5.2f) = %f\n",n,X[x],Interpolation(n,X[x]));
        printf("\nenter the exponential you want:");

    }
    return 0;
}
