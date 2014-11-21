#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CS.h"



double DC(int k, double W[nnest][nobj]){
    
    double sumTer1 = 0;
    double sumTer3 = 0;
    int i;
    for (i = 0; i < nnest; i++)
    {
        double prodTer1 = 1;
        int j;
        for (j = 0; j < k; j++)
        {
            prodTer1 = prodTer1 * (1 + ((1.0/2.0) * fabs(W[i][j] - 1.0/2.0)) - ((1.0/2.0) * pow(fabs(W[i][j] - 1.0/2.0), 2)) );
        }
        sumTer1 = sumTer1 + prodTer1;
        double sumTer2 = 0;
        int r;
        for (r = 0; r < nnest; r++)
        {
            double prodTer2 = 1;
            int j;
            for (j= 0; j < k; j++)
            {
                prodTer2 = prodTer2 * (1 + ((1.0/2.0) * fabs(W[i][j] - 1.0/2.0)) + ((1.0/2.0) * fabs(W[r][j] - 1.0/2.0)) - ((1.0/2.0) * fabs(W[i][j] - W[r][j])) );
            }
            sumTer2 = sumTer2 + prodTer2;
        }
        sumTer3 = sumTer3 + sumTer2;
    }

    double DC = sqrt((pow(13.0/12.0,k) - (2.0/nnest * sumTer1)) + (1.0/pow(nnest,2)) * sumTer3);
    
    return DC;

}

void generarArchivo(int k, double W[nnest][nobj]){
    /* Se genera un archivo .dem para leerlo con el comando load "pesos.dem" en gnuplot para modelar los puntos*/
    FILE *f;
    f = fopen("pesos.dem", "w");
    fprintf(f, "splot '-'\n");
    int i;
    for (i = 0; i < nnest; i++)
    {
        int j;
        for (j = 0; j < k; j++)
        {
            fprintf(f,"%f\t",W[i][j]);
        }
        fprintf(f, "\n" );
    }
    fprintf(f, "e" );
    fclose(f);
    printf("\n Se genero el archivo correctamente\n");
}

void primos(int k, int num, double primo[nobj])
{
    int n, p, i;

    n = 2;
    i = 1;
    while(i <= k)
    {
        for(p = 2; n % p != 0; p++);
            if(p == n)
            {
                primo[i-1] = n;
                i++;
            }
        n++;
    }
} 

void weightVectors(int k, double W[nnest][nobj], double primo[nobj]){
   double U[nnest][k];
    double n2veces = 2*nnest;

    int i;
    for (i = 0; i < nnest; i++)
    {

        U[i][0] = (2*(i+1)-1)/n2veces;
        int j;
        for (j = 1; j < k-1; j++)
        {
            U[i][j] = 0;
            double f = (double)1.0/primo[j-1];
            double D = i+1.0;
            while(D>0){
                U[i][j] = U[i][j] + f * fmod(D, primo[j-1]);
                D = floor(D / primo[j-1]);
                f = f/primo[j-1];
            }
    
        }   
        
    }

    int t;
    for (t = 0; t < nnest; t++){
        double divTemp;
        int j;
        for (i = 0; i < k-1; i++){
            
            double prodTemp = 1;
            
                for (j = 0; j < (i+1)-1; j++)
                {
                    divTemp = k-(j+1);
                    prodTemp = prodTemp * pow(U[t][j], 1.0/divTemp );

                }
            divTemp = k-(i+1);
            W[t][i] = (1 - pow( U[t][i], (1.0/divTemp) ) ) * (prodTemp);
            
        }

        double prodTemp = 1;
        for (j = 0; j < (i+1)-1; j++)
        {
            divTemp = k-(j+1);
            prodTemp = prodTemp * pow(U[t][j], 1.0/divTemp );

        }
        W[t][k-1] = prodTemp;
    } 
}

int main(){
    int k = nobj;
    double primo[k], W[nnest][nobj];
    primos(k,nnest, primo);
    weightVectors(k,W,primo);
    
    int i;
    for (i = 0; i < nnest; i++)
    {
        int j;
        
        for (j = 0; j < k; j++)
        {
            printf("%f\t",W[i][j]);
        }
        printf("\n");
    }
    printf("Medida de uniformidad Discrepancia-L2 Centrada: %f",DC(k, W));
    generarArchivo(k,W);
    return 0;
}
