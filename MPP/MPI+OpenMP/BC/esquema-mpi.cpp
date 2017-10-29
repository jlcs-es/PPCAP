#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <mpi.h>
#include <signal.h>
#include <unistd.h>

//Spanish Parallel programming Contest 2017. Problem B.
//Multiplication of tailed matrices. Version MPI.
//Schema for In/Out, validation and execution time

void generar(double *m, int t) {
    int i;
    for (i = 0; i < t; i++) {
        m[i] = (20. * rand()) / RAND_MAX-10.;
    }
}

void generatezerosrows(int n,double *a) {
    int pos;
    for(int i=0;i<n;i++)
    {
        pos=n*((1.*rand())/RAND_MAX);
        for(int j=pos;j<n;j++)
            a[i*n+j]=0.;
    }
}

void generatezeroscolumns(int n,double *a) {
    int pos;
    for(int i=0;i<n;i++) {
        pos=n*((1.*rand())/RAND_MAX);
        for(int j=pos;j<n;j++)
            a[j*n+i]=0.;
    }
}

void escribir(double *m, int t) {
    int i, j;
    for (i = 0; i < t; i++) {
        for (j = 0; j < t; j++)
            printf("%.4lf ", m[i * t + j]);
        printf("\n");
    }
    printf("\n");
}

void escribirresult(double *a,int N,int salida) {
    int i;
    for(i=0;i<N;i++) {
        if((i%salida)==0) {
            printf("%lf \n",a[i]);
        }
    }
    printf("\n");
}

/*
 mseconds - returns elapsed milliseconds since Jan 1st, 1970.
*/
long long mseconds() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec*1000 + t.tv_usec/1000;
}

static void alarm_handler(int sig) {
    fprintf(stderr, "Time Limit Exceeded\n");
    abort();
}

extern void sec(int t,double *a,double *b,double *c,int nodo,int np);

int main(int argc,char *argv[]) {
    int i,N;
    int cuantos,semilla,salida;
    long long ti,tf,tt=0;
    double *a,*b,*c;
    int nodo,np;

    FILE *stats_file = fopen("stats", "w");
    struct sigaction sact;
    sigemptyset(&sact.sa_mask);
    sact.sa_flags = 0;
    sact.sa_handler = alarm_handler;
    sigaction(SIGALRM, &sact, NULL);
    alarm(40); /* time limit */

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&nodo);

    // The number of test cases is read
    // and sent to all the processes
    if(nodo==0) {
        scanf("%d",&cuantos);
        MPI_Bcast(&cuantos,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    else {
        MPI_Bcast(&cuantos,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    for(i=0;i<cuantos;i++) {
        if(nodo==0) {
            scanf("%d",&N); // Matrices size
            scanf("%d",&semilla); // seed for random generation
            scanf("%d",&salida); // to determine the elements to be written to the output
            // Space for the matrix, the values, rows and columns
            a = (double *) calloc(sizeof(double),N*N);
            b = (double *) calloc(sizeof(double),N*N);
            c = (double *) calloc(sizeof(double),N*N);
            srand(semilla);
            generar(a,N*N);
            generatezerosrows(N,a);
            generar(b,N*N);
            generatezeroscolumns(N,b);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        ti=mseconds();

        sec(N,a,b,c,nodo,np);

        MPI_Barrier(MPI_COMM_WORLD);
        
        tf=mseconds();
        if(nodo==0) {
            if(i!=0) tt+=tf-ti;
            escribirresult(c,N*N,salida);
            free(a);
            free(b);
            free(c);
        }
    }
    if(nodo==0) {
        fprintf(stats_file, "%Ld\n", tt);
        fclose(stats_file);
    }
    MPI_Finalize();
    return 0;
}

