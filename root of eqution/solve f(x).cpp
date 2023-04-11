#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;

double solve (double x){
    //return exp(x) - 5.6*x*sin(x) + 3.2*x + 1.8; /*第一題*/
    return  1 + pow(x,2) - tan(x); /*第二題*/
}

double dsolve (double x){
    //return exp(x) -5.6*x*cos(x) -5.6*sin(x) + 3.2; /*第一題*/
    return  2*x - 1/ (cos(x)*cos(x)) ; /*第二題*/
}

double gxsolve (double x){
    //return (exp(x) + 3.2*x + 1.8) / (5.6*sin(x)); /*第一題*/
    return  sqrt(tan(x)-1) ; /*第二題*/
}

void Bisection (double l,double r){
    while(r-l >= 0.000002){
        double m = (l+r)/2;
        if(solve(l)*solve(m)<0)
            r = m;
        else
            l = m;
        printf("root is between %f and %f\n",l,r);
    }
    printf("root is  %f \n",(l+r)/2);
    printf("\n");
}

void FalsePosition (double a,double b){
    double old = 0, m = 1;

    while(abs(m-old) >= 0.000001){
        old = m;
        m = ( a*solve(b) - b*solve(a) ) / ( solve(b)- solve(a) );
        if(solve(a)*solve(m)<0)
            b = m;
        else
            a = m;
        printf("root is between %f and %f\n",a,b);
    }
    printf("root is  %f \n",m);
    printf("\n");
}

void ModifyFalsePosition (double a,double b){
    double old = 0, m = 1, Fa = solve(a), Fb = solve(b),Fm;

    while(abs(m-old) >= 0.000001){
        old = m;
        m = ( a*Fb - b*Fa ) / ( Fb-Fa );
        Fm = solve(m);
        if(Fa*Fm < 0){
            b = m;
            Fb = solve(b);
            Fa = Fa/2;
        }

        else{
            a = m;
            Fa = solve(a);
            Fb = Fb/2;
        }
        printf("root is between %f and %f\n",a,b);
    }
    printf("root is  %f \n",m);
    printf("\n");
}

void NewtonMethod (double a){
    double delta = -1*solve(a)/dsolve(a);
    int cnt = 0;

    while(abs(delta) >= 0.000001){
        a = a + delta;
        printf("root is approaching at %f\n",a);
        delta = -1*solve(a)/dsolve(a);
        cnt++;
        if(cnt > 40){
            printf("need to change the guess better!!\n");
            break;
        }
    }
    printf("root is  %f \n",a);
    printf("\n");
}

void SecantMethod (double a,double b){
    double x = ( a*solve(b) - b*solve(a) ) / ( solve(b)- solve(a) );
    int cnt = 0;

    while(abs(b-a) >= 0.000001){
        a = b;
        b = x;
        printf("root is approaching at %f\n",b);
        x = ( a*solve(b) - b*solve(a) ) / ( solve(b)- solve(a) );
        cnt++;
        if(cnt > 40){
            printf("need to change the guess better!!\n");
            break;
        }
    }
    printf("root is  %f \n",b);
    printf("\n");
}

void FixedPointMethod (double x0){
    double x = gxsolve(x0);
    int cnt = 0;

    while(abs(x0-x) > 0.000001){
        printf("%f %f\n",x0,x);
        x0 = x;
        x = gxsolve(x0);
        printf("root is approaching at %f\n",x);
        cnt++;
        if(cnt > 40){
            printf("need to change the guess better!!\n");
            break;
        }
    }
    printf("root is  %f \n",x);
    printf("\n");
}

int main()
{   for(double i = -30; i <30; i+=1){
        if(solve(i)*solve(i+1) < 0){
            printf("root is between %f and %f\n",i,i+1);
            printf("Using Bisection Method:\n");
            Bisection(i,i+1);

            printf("Using False Position:\n");
            FalsePosition(i,i+1);

            printf("Using Modify False Position:\n");
            ModifyFalsePosition(i,i+1);
        }
    }

    double guess = -15;
    printf("Using Newton Method and guess x is -15:\n");
    NewtonMethod(guess);
    guess = 4;
    printf("Using Newton Method and guess x is 4:\n");
    NewtonMethod(guess);

    double guessleft = -20, guessright = -15;
    printf("Using Secant Method and guess x is between -20 and -15:\n");
    SecantMethod(guessleft,guessright);
    guessleft = 1, guessright = 8;
    printf("Using Secant Method and guess x is between 1 and 8:\n");
    SecantMethod(guessleft,guessright);

    double guessx = -2;
    printf("Using Fixed Point Method and guess x is -2:\n");
    FixedPointMethod(guessx);

    return 0;
}
