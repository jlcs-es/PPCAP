#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>


int NUM_THREADS = 1;
 

void mms(double *a, int fa, int ca, int lda, double *b, int fb, int cb, int ldb, double *c, int fc, int cc, int ldc) {
    int i, j, k;
    double s;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for private(i, j, k, s) schedule(dynamic)
    for (i = 0; i < fa; i++) 
        for (j = 0; j < cb; j++) {
            s = 0.;
            for (k = 0; k < ca; k++)
                s += a[i * lda + k] * b[k * ldb + j];
            c[i * ldc + j] = s;
        }
}


void matriz_matriz_mpi(int n,double *a, double *b,double *c,int nodo, int np) {
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    int i,fa=n,ca=n,lda=n,fb=n,cb=n,ldb=n,fc=n,cc=n,ldc=n,position, prevposition;
    double s;

    if(nodo!=0)
    {
        a=(double *) malloc(sizeof(double)*n*n/np);
        b=(double *) malloc(sizeof(double)*n*n);
        c=(double *) malloc(sizeof(double)*n*n/np);
    }

    char * buffer = (char *) malloc(sizeof(double)*(n*n + n*n/np));

    if (nodo == 0) {
        position = 0;
        MPI_Pack(b, fb*cb, MPI_DOUBLE, buffer, sizeof(double)*(n*n + n*n/np), &position, MPI_COMM_WORLD);
        prevposition = position;
        for (i = 1; i < np; i++){
            MPI_Pack(&a[i * lda * fa / np], fa / np * ca, MPI_DOUBLE, buffer, sizeof(double)*(n*n + n*n/np), &position, MPI_COMM_WORLD);
            MPI_Send(buffer, sizeof(double)*(n*n + n*n/np), MPI_PACKED, i, 20, MPI_COMM_WORLD);
            position = prevposition;
        }
    } else {
        MPI_Recv(buffer, sizeof(double)*(n*n + n*n/np), MPI_PACKED, 0, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        position = 0;
        MPI_Unpack(buffer, sizeof(double)*(n*n + n*n/np), &position, b, fb*cb, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack(buffer, sizeof(double)*(n*n + n*n/np), &position, a, fa / np * ca, MPI_DOUBLE, MPI_COMM_WORLD);
    }
    free(buffer);

    mms(a, fa / np, ca, lda, b, fb, cb, ldb, c, fc / np, cc, ldc);
    MPI_Gather(c, fc / np * cc, MPI_DOUBLE,c, fc / np * cc, MPI_DOUBLE,0, MPI_COMM_WORLD);
   if(nodo!=0)
   {
       free(a);
       free(b);
       free(c);
   }
}

void generar_vector(double *v,int n,double l,double u)
{
  int i;

  for(i=0;i<n;i++)
    v[i]=((double) rand()/RAND_MAX)*(u-l)+l;
}

int main(int argc,char *argv[])
{

  double *a,*b,*c,l,u;
  int nodo,np,i,j,N; 

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&nodo);


  int ti1,tf1,si,sf;
  struct timeval *tv;
  struct timezone *tz;
  tv=(struct timeval *) malloc(sizeof(struct timeval));
  tz=(struct timezone *) malloc(sizeof(struct timezone));

  if(nodo==0)
    {
            printf("De el tamaño N :");
            scanf("%d",&N);
            printf("De el numero de hilos openmp: ");
            scanf("%d",&NUM_THREADS);
            printf("De los valores inferior y superior: ");
            scanf("%lf %lf",&l,&u);

            a=(double *) malloc(sizeof(double)*N*N);
            b=(double *) malloc(sizeof(double)*N*N);
            c=(double *) malloc(sizeof(double)*N*N);

            generar_vector(a,N*N,l,u);
            generar_vector(b,N*N,l,u);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    gettimeofday(tv,tz);
    si=(tv->tv_sec);
    ti1=(tv->tv_usec);
    ti1=si*1000000+ti1;

    matriz_matriz_mpi(N, a, b, c, nodo, np);

    MPI_Barrier(MPI_COMM_WORLD);

    gettimeofday(tv,tz);
    sf=(tv->tv_sec);

    tf1=(tv->tv_usec);
    tf1=sf*1000000+tf1;
    if(nodo == 0){
        printf("Tamaño matriz %d - MPI NP %d - OMP NT %d : %.6f seg\n",N,np,NUM_THREADS,(tf1-ti1)/1000000.);
        free(a);
        free(b);
        free(c);
    }
    
    free(tv);
    free(tz);
    
  MPI_Finalize();
  return 0;
}

