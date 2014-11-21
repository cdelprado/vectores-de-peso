#define iter 1000

#define nnest 120
#define nobj 3
#define d 15

typedef struct Nest{
    double v[d];
    double fx[2];
}Nests;

void CS(int, int, int, Nests *,double  *(*fobj)(double[],int));

double *(*fobj)(double x[], int);
