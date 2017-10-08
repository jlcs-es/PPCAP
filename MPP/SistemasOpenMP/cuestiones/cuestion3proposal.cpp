// TODO

//En linux, 
//compilar con: g++ -O3 cuestion3proposal.cpp -fpermissive -lpthread
//y la -D DEBUG si se quiere compilar la zona de escritura de datos: g++ -O3 ejemplofork.cpp -D DEBUG
//ejecutar con: ./a.out


#include <pthread.h>

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


#define NUM_THREADS 8

struct thread_data{
    int arows;
    int acols_brows;
    int bcolumns;
    double *a;
    double *b;
    double *c;
};


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

void addMatrices(double *a, double *b, int t, double *c) {
    for(int i = 0; i < t; i++)
        c[i] = a[i]+b[i];
}

void *multiply(struct thread_data *data) {
    double *a       = data->a;
    double *b       = data->b;
    int arows       = data->arows;
    int acols_brows = data->acols_brows;
    int bcolumns    = data->bcolumns;
    double *c       = data->c;
    multiplyMatrices(a, b, arows, acols_brows, bcolumns, c);
}

int main() {
    int N;
    int semilla;
    long long ti,tf;
    double *a,*b,*c, *csec;
    void *status; //se usa para comprobar el estado de los hilos
    pthread_t threads[NUM_THREADS]; //se declara un array de hilos para ponerlos en marcha
    struct thread_data thread_data_array[NUM_THREADS];

    //cin >> N; //se lee el tamaño de la matriz
    N = 1000; // 
    cin >> semilla; // se lee la semilla
    srand(semilla);

    a = (double *) calloc(sizeof(double),N*N);
    b = (double *) calloc(sizeof(double),N*N);
    c = (double *) calloc(sizeof(double),N*N);
    csec = (double *) calloc(sizeof(double),N*N);

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

    // Multiplicación paralela
    double *A11 = (double *) calloc(sizeof(double),N*N/4.);
    double *A12 = (double *) calloc(sizeof(double),N*N/4.);
    double *A21 = (double *) calloc(sizeof(double),N*N/4.);
    double *A22 = (double *) calloc(sizeof(double),N*N/4.);
    double *B11 = (double *) calloc(sizeof(double),N*N/4.);
    double *B12 = (double *) calloc(sizeof(double),N*N/4.);
    double *B21 = (double *) calloc(sizeof(double),N*N/4.);
    double *B22 = (double *) calloc(sizeof(double),N*N/4.);
    copiar(a, A11, N, 0, 0, N/2, N/2);
    copiar(a, A12, N, 0, N/2, N/2, N/2);
    copiar(a, A21, N, N/2, 0, N/2, N/2);
    copiar(a, A22, N, N/2, N/2, N/2, N/2);
    copiar(b, B11, N, 0, 0, N/2, N/2);
    copiar(b, B12, N, 0, N/2, N/2, N/2);
    copiar(b, B21, N, N/2, 0, N/2, N/2);
    copiar(b, B22, N, N/2, N/2, N/2, N/2);
    
    
    
    ti=mseconds();

    //en este bucle se crearán los hilos
    for(int i=0; i<NUM_THREADS; i++){
        thread_data_array[i].arows          = N/2;
        thread_data_array[i].acols_brows    = N/2;
        thread_data_array[i].bcolumns       = N/2;
        thread_data_array[i].a              = &a[i*(N/NUM_THREADS)*N];
        thread_data_array[i].b              = b;
        thread_data_array[i].c              = &c[i*(N/NUM_THREADS)*N];
        //salvo el último, que ordena todos los datos restantes
        if(i==(NUM_THREADS-1))
        {
            thread_data_array[i].arows      = N-i*(N/NUM_THREADS);
        }

        //se crea el hilo pasando su dirección, la función que ejecuta, y los parámetros de la función en una estructura
        pthread_create(&threads[i], NULL, multiply,(void *) &thread_data_array[i]);
    }

    //el hilo maestro espera a que los esclavos acaben
    for(int i=0; i<NUM_THREADS; i++) {
        pthread_join(threads[i], &status);
    }

    tf=mseconds();

    #ifdef DEBUG
        escribir(c, N, N);
    #endif
    
    cout <<"Tiempo paralelo: " <<tf-ti<<" microsegundos"<<endl;
    cout <<"Diferencia: "<< comparar(csec, c, N*N) <<endl;

    pthread_exit(NULL);

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



void copiar(double* a, double *b, int n_cols, int row_init, int column_init, int row_end, int column_end) {
    for(int i = row_init; i<row_end; i++) {
        for(int j = column_init; j < column_end; j++) {
            b[ (i-row_init)*(row_end-row_init)+(j-column_init) ] = a[i*n_cols + j];
        }
    }
}