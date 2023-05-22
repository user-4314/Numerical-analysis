#include <iostream>
#include <cstdio>
using namespace std;
const double PI = 3.1415926535897932;
double solvef(double x){
    return 1/(x*x+1);
}

void Trapezoidal(){
    double left = 0,right = 1;

    for(int n = 10;n <= 10000000;n*=10){
        double delta = ( right - left ) / n, sum = 0;
        sum = sum + solvef(left) + solvef(left + n*delta);

        for(int i = 1;i < n;i++)
            sum = sum + 2 * solvef(left + i*delta);

        sum = sum * delta / 2;
        printf("when n = %8d ,using Trapezoidal Rule answer is %.16f\n",n,sum*4);
    }
}

void SimpsonOneThird(){
    double left = 0,right = 1;

    for(int n = 10;n <= 10000;n*=10){
        double delta = ( right - left ) / n, sum = 0;
        sum = sum + solvef(left) + solvef(left + n*delta);

        for(int i = 1;i < n;i++){
            if(i%2 == 0)
                sum = sum + 2 * solvef(left + i*delta);
            else
                sum = sum + 4 * solvef(left + i*delta);
        }
        sum = sum * delta / 3;
        printf("when n = %8d ,using Simpson 1/3 Rule answer is %.16f\n",n,sum*4);
    }
}

void SimpsonThreeEighths(){
    double left = 0,right = 1;

    for(int n = 3;n <= 300;n*=10){
        double delta = ( right - left ) / n, sum = 0;
        sum = sum + solvef(left) + solvef(left + n*delta);

        for(int i = 1;i < n;i++){
            if(i%3 == 0)
                sum = sum + 2 * solvef(left + i*delta);
            else
                sum = sum + 3 * solvef(left + i*delta);
        }
        sum = sum * delta * 3 / 8;
        printf("when n = %8d ,using Simpson 3/8 Rule answer is %.16f\n",n,sum*4);
    }
}

int main()
{
    printf("using Trapezoidal Rule\n");
    Trapezoidal();
    printf("                                      the real PI = %.16f\n",PI);

    printf("using Simpson 1/3 Rule\n");
    SimpsonOneThird();
    printf("                                      the real PI = %.16f\n",PI);

    printf("using Simpson 3/8 Rule\n");
    SimpsonThreeEighths();
    printf("                                      the real PI = %.16f\n",PI);

    return 0;
}
