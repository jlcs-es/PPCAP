#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <omp.h>
#include <signal.h>
#include <unistd.h>

// Ejemplo de multiplicación de matrices
//esquema para la entrada/salida, validación y toma de tiempos
//para el proyecto de Iniciación al paralelismo en módulos de FP
//Universidad de Murcia-Instituto de Enseñanza Secundaria Ingeniero de la Cierva, Murcia, 2014

void generar(double *m, int t) {
  int i;

  for (i = 0; i < t; i++) {
      m[i] = (20. * rand()) / RAND_MAX-10.;
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

void escribirresult(double *a,int N,int salida)
{
  int i;

  for(i=0;i<N;i++)
  {
        if((i%salida)==0)
        {
          printf("%lf \n",a[i]);
        }
  }
  printf("\n");
}

/*
c
c     mseconds - returns elapsed milliseconds since Jan 1st, 1970.
c
*/
long long mseconds(){
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec*1000 + t.tv_usec/1000;
}

static void alarm_handler(int sig) {
  fprintf(stderr, "Time Limit Exceeded\n");
  abort();
}

extern void sec(int t,double *a,double *b,double *c);

int main(int argc,char *argv[]) {
  int i,N;
  int cuantos,semilla,salida;
  long long ti,tf,tt=0;
  double *a,*b,*c;

  FILE *stats_file = fopen("stats", "w");

  struct sigaction sact;
  sigemptyset(&sact.sa_mask);
  sact.sa_flags = 0;
  sact.sa_handler = alarm_handler;
  sigaction(SIGALRM, &sact, NULL);
  alarm(60);  /* time limit */
 
  scanf("%d",&cuantos);

  for(i=0;i<cuantos;i++)
  {
      scanf("%d",&N);                   // Primero se da el tamaño de las matrices
      scanf("%d",&semilla);             // después la semilla
      scanf("%d",&salida);                 // y un valor para seleccionar qué elementos se escriben

// Space for the matrix, the values, rows and columns
      a = (double *) calloc(sizeof(double),N*N);
      b = (double *) calloc(sizeof(double),N*N);
      c = (double *) calloc(sizeof(double),N*N);

      srand(semilla);
      generar(a,N*N);
      generar(b,N*N);

    ti=mseconds(); 
    sec(N,a,b,c);
    tf=mseconds(); 
      if(i!=0) tt+=tf-ti;
      escribirresult(c,N*N,salida);
      free(a);
      free(b);
      free(c);
  }
  
    fprintf(stats_file, "%Ld\n", tt);
    fclose(stats_file);
  return 0;
}
