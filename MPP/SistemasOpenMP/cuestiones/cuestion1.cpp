//En linux, 
//compilar con: g++ -O3 cuestion1.cpp
//y la -D DEBUG si se quiere compilar la zona de escritura de datos: g++ -O3 ejemplofork.cpp -D DEBUG
//ejecutar con: ./a.out

#include <stdio.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include <time.h>
using namespace std;




void generar(double *m, int t);
void escribir(double *m, int r, int c);
void escribirresult(double *a,int N,int salida);
double comparar(double *x,double *y,int t);
long long mseconds();


void multiplyMatrices(double *a, double *b, int arows, int acols_brows, int bcolumns, double *c) {
    double s;
    for (int i = 0; i < arows; i++)
        for (int j = 0; j < bcolumns; j++) {
            s = 0.;
            for (int k = 0; k < acols_brows; k++)
                s += a[i * acols_brows + k] * b[k * bcolumns + j];
            c[i * bcolumns + j] = s;    // c dim = arows * bcolumns
    }
}

int main() {
    int N;
    int semilla;
    long long ti,tf;
    int status;
    double *a,*b,*c, *csec;

    pid_t pid,cpid; //identificador de proceso

    cin >> N; //se lee el tamaño de la matriz
    cin >> semilla; // se lee la semilla
    srand(semilla);

    a = (double *) calloc(sizeof(double),N*N);
    b = (double *) calloc(sizeof(double),N*N);
    csec = (double *) calloc(sizeof(double),N*N);
    c = (double *) mmap(NULL, sizeof(double)*N*N, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);

    generar(a, N*N);
    generar(b, N*N);

    // Multiplicación secuencial
    ti=mseconds(); 
    multiplyMatrices(a, b, N, N, N, csec);
    tf=mseconds();
    
    cout <<"Tiempo secuencial: " <<tf-ti<<" microsegundos"<<endl;

    //Si se compila con -D DEBUG se escriben los datos
    #ifdef DEBUG
        escribir(csec, N, N);
    #endif

    //toma de tiempos de la multiplicación paralela
    ti=mseconds();

    //fork devuelve 0 al proceso hijo
    //y el número de proceso creado al padre
    pid = fork();
 
    //El padre hace esta parte
    if (pid != 0) {
        multiplyMatrices(a, b, N/2, N, N, c);
        //El padre espera a que el hijo acabe
        while((cpid=wait(&status))>0);
    } 
    //y el hijo hace esta parte
    else {
        multiplyMatrices(&a[N*N/2], b, N/2, N, N, &c[N*N/2]);
        exit(0);
    }
    tf=mseconds();
    #ifdef DEBUG
        escribir(c, N, N);
    #endif

    cout <<"Tiempo paralelo: " <<tf-ti<<" microsegundos"<<endl;
    cout <<"Diferencia: "<< comparar(csec, c, N*N) <<endl;
}




long long mseconds() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec*1000 + t.tv_usec/1000;
}

void generar(double *m, int t) {
    int i;

    for (i = 0; i < t; i++) {
        m[i] = (20. * rand()) / RAND_MAX-10.;
    }
}
  
void escribir(double *m, int r, int c) {
    int i, j;

    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++)
        printf("%.4lf ", m[i * c + j]);
        printf("\n");
    }
    printf("\n");
}


double comparar(double *x,double *y, int t) {
    int i;
    double total=0;

    for(i=0;i<t;i++)
    {
        total=(y[i]>x[i]?y[i]-x[i]:x[i]-y[i]);
    }
    return total;
}
