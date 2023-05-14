#include <iostream>
#include <cstdio>
using namespace std;
const int xN = 9,yN = 8;
double X[xN] ={0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7},Y[yN] = {0.25,0.37,0.49,0.61,0.73,0.85,0.97,1.09};
double Z[xN][yN] =
{
6.9975, 7.7967 , 8.6823 , 9.6543 , 10.7127 , 11.8575 , 13.0887 , 14.4063 ,
6.1775, 6.9287 , 7.7663 , 8.6903 ,  9.7007 , 10.7975 , 11.9807 , 13.2503 ,
5.4375, 6.1407 , 6.9303 , 7.8063 ,  8.7687 ,  9.8175 , 10.9527 , 12.1743 ,
4.7775, 5.4327 , 6.1743 , 7.0023 ,  7.9167 ,  8.9175 , 10.0047 , 11.1783 ,
4.1975, 4.8047 , 5.4983 , 6.2783 ,  7.1447 ,  8.0975 ,  9.1367 , 10.2623 ,
3.6975, 4.2567 , 4.9023 , 5.6343 ,  6.4527 ,  7.3575 ,  8.3487 ,  9.4263 ,
3.2775, 3.7887 , 4.3863 , 5.0703 ,  5.8407 ,  6.6975 ,  7.6407 ,  8.6703 ,
2.9375, 3.4007 , 3.9503 , 4.5863 ,  5.3087 ,  6.1175 ,  7.0127 ,  7.9943 ,
2.6775, 3.0927 , 3.5943 , 4.1823 ,  4.8567 ,  5.6175 ,  6.4647 ,  7.3983
};
double xFDT[xN][xN] = {0}, yFDT[xN][yN][yN] = {0};

void MakingXTable(double y){
    for(int i=1;i<xN;i++){
        for(int j=0;j<xN-i;j++)
            xFDT[i][j]= xFDT[i-1][j+1]-xFDT[i-1][j];
    }
    printf("F(x,%.2f)different table:\n",y);
    for(int i=0;i<xN;i++){
        for(int j=0;j<xN-i;j++)
            printf("%10f ",xFDT[i][j]);
        printf("\n");
    }
    printf("\n");
}

void MakingYTable(int n){
    for(int i=0;i<yN;i++)
        yFDT[n][0][i] = Z[n][i];
    for(int i=1;i<yN;i++){
        for(int j=0;j<yN-i;j++)
            yFDT[n][i][j]= yFDT[n][i-1][j+1]-yFDT[n][i-1][j];
    }
    printf("F(%.2f,y)different table:\n",X[n]);
    for(int i=0;i<yN;i++){
        for(int j=0;j<yN-i;j++)
            printf("%10f ",yFDT[n][i][j]);
        printf("\n");
    }
    printf("\n");
}

void yForwardDifferenceTable(double x,int n){
    double s = (x-Y[0]) / abs( (Y[1]-Y[0]) );
    double y = yFDT[n][0][0];
    for(int i=1;i<yN;i++){
        double P = 1;
        for(int j=0;j<i;j++)
            P = P * (s-j) / (j+1);
        y = y + P*yFDT[n][i][0];
    }
    xFDT[0][n] = y;
    printf("bring F(%.2f,%.2f) by Forward Difference Table then z = %f\n",X[n],x,y);
}

void yBackwardDifferenceTable(double x,int n){
    double s = (x-Y[yN-1]) / abs( (Y[1]-Y[0]) );
    double y = yFDT[n][0][yN-1];
    for(int i=1;i<yN;i++){
        double P = 1;
        for(int j=0;j<i;j++)
            P = P * (s+j) / (j+1);
        y = y + P*yFDT[n][i][yN-i-1];
    }
    xFDT[0][n] = y;
    printf("bring F(%.2f,%.2f) by Backward Difference Table then z = %f\n",X[n],x,y);
}

void yNewtonForward(double x,int n){
    int pos = yN/2;
    double s = (x-Y[pos]) / abs( (Y[1]-Y[0]) ) -1;
    double y = yFDT[n][0][pos];
    for(int i=1;pos>0 && i<yN;i++){
        double P = 1;
        if(i%2 ==0)
            pos--;
        else
            s++;
        for(int j=0;j<i;j++)
            P = P * (s-j) / (j+1);
        y = y + P*yFDT[n][i][pos];
    }
    xFDT[0][n] = y;
    printf("bring F(%.2f,%.2f) by Newton Forward then z = %f\n",X[n],x,y);
}

void xForwardDifferenceTable(double x,double y){
    double s = (x-X[0]) / abs( (X[1]-X[0]) );
    double z = xFDT[0][0];
    for(int i=1;i<xN;i++){
        double P = 1;
        for(int j=0;j<i;j++)
            P = P * (s-j) / (j+1);
        z = z + P*xFDT[i][0];
    }
    printf("bring F(%.2f,%.2f) by Forward Difference Table then z = %f\n",x,y,z);
}

int main()
{
    for(int i=0;i<xN;i++)
        MakingYTable(i);

    printf("1. solve f(0.6,0.7):\n");
    printf("\n");

    for(int i=0;i<xN;i++)
        yNewtonForward(0.7,i);
    printf("\n");
    MakingXTable(0.7);
    xForwardDifferenceTable(0.6,0.7);
    printf("\n");


    printf("2. solve f(0.24,0.9):\n");
    printf("\n");

    for(int i=0;i<xN;i++)
        yBackwardDifferenceTable(0.9,i);
    printf("\n");
    MakingXTable(0.9);
    xForwardDifferenceTable(0.24,0.9);
    printf("\n");


    printf("3. solve f(0.42,0.48):\n");
    printf("\n");

    for(int i=0;i<xN;i++)
        yForwardDifferenceTable(0.48,i);
    printf("\n");
    MakingXTable(0.48);
    xForwardDifferenceTable(0.42,0.48);
    printf("\n");
    return 0;
}
