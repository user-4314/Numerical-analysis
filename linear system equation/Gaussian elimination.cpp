#include <iostream>
#include <cstdio>
using namespace std;
const int N = 25;
float a[N][N+1],b[N][N+1];
float x[N] = {0},xb[N] = {0};

void PartialPivoting(float a[N][N+1],int pos){
    float ma = a[pos][pos];
    int mapos = pos;
    for(int i = pos+1;i < N;i++){
        if(abs(a[i][0]) > ma){
            ma = abs(a[i][0]);
            mapos = i;
        }
    }
    for(int j = 0;j < N+1;j++){
        float counter = a[mapos][j];
        a[mapos][j] = a[pos][j];
        a[pos][j] = counter;
    }
}

void PrintArray(float a[N][N+1]){
    for(int i = 0;i < N;i++){
        for(int j = 0;j < N+1;j++)
            printf("%8.2f ",a[i][j]);
        printf("\n");
    }
    printf("\n");
}

int main()
{
    for(int i = 0;i < N;i++){
        for(int j = 0;j < N;j++){
            a[i][j] = rand()%30 -10;
            b[i][j] = a[i][j];
        }
    }
    for(int i = 0;i < N;i++){
        a[i][N] = rand()%55 -5;
        b[i][N] = a[i][N];
    }
    printf("random produce linear system equation:\n");
    PrintArray(a);

    for(int k = 0;k < N-1;k++){
        PartialPivoting(a,k);
        for(int i = k+1;i < N;i++){
            float d = a[i][k]/a[k][k];
            for(int j = k;j <= N;j++){
                a[i][j] = a[i][j] - d*a[k][j];
            }
        }
    }
    printf("after Gaussian elimination:\n");
    PrintArray(a);

    for(int k = N-1;k>=0;k--){
        float sum = 0;
        for(int j = k+1;j < N;j++){
            sum = sum + a[k][j]*x[j];
        }
        x[k] = (a[k][N] - sum) / a[k][k];
    }
    printf("the solution:\n");
    for(int i = 0;i < N;i++)
        printf("x%d is %f\n",i,x[i]);
    printf("\n");

    /*printf("check the answer correct or not\n");
    for(int i=0;i<N;i++){
        float check = 0;
        for(int j=0;j<N;j++)
           check = check + b[i][j] * x[j];
        printf("in line %d %f = %f\n",i,check,b[i][N]);
    }
    printf("\n");*/

    for(int k = 0;k < N;k++){
        for(int i = 0;i < N;i++){
                float d;
                if(i != k)
                     d = b[i][k]/b[k][k];
                else
                    continue;

                for(int j = k;j <= N;j++)
                    b[i][j] = b[i][j] - d*b[k][j];
        }
    }
    printf("after Gauss Jordon method:\n");
    PrintArray(b);

    for(int k = N-1;k>=0;k--)
        xb[k] = b[k][N] / b[k][k];
    printf("the solution:\n");
    for(int i = 0;i < N;i++)
        printf("x%d is %f\n",i,xb[i]);
    printf("\n");
    return 0;
}
