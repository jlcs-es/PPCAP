#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <signal.h>
#include <unistd.h>

//Spanish Parallel programming Contest 2017. Problem I - heterosolar 
//Maximum coincidence with a mask in 2D. CUDA version.
//Schema for In/Out, validation and execution time

void generar(char *m, int t,int sup) {
  int i;

  for (i = 0; i < t; i++) {
      m[i] = (char) (((1. * rand()) / RAND_MAX)*sup)+'a';
  }
}

void escribir(char *m, int t) {
  int i;

  for (i = 0; i < t; i++) {
      printf("%c ", m[i]);
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

extern int sec(int,char *,int,char *);

int main(int argc,char *argv[]) {
  int N,M,cuantos;
  bool correcto=true;
  int semilla,upper;
  char *A,*B;
  long long ti,tf,tt=0;

  //FILE *stats_file = fopen("stats", "w");

  struct sigaction sact;
  sigemptyset(&sact.sa_mask);
  sact.sa_flags = 0;
  sact.sa_handler = alarm_handler;
  sigaction(SIGALRM, &sact, NULL);
  alarm(40);  /* time limit */

  scanf("%d",&cuantos);

  for(int i=0;i < cuantos;i++)
  {
      scanf("%d",&N);                   // Matrices size
      scanf("%d",&M);                   // mask size
      scanf("%d",&semilla);             // seed for random generation
      scanf("%d",&upper);                 // upper value for random generation

// Space for the matrix, the values, rows and columns
      A = (char *) malloc(sizeof(double)*N*N);
      B = (char *) malloc(sizeof(double)*M*M);

      srand(semilla);
     
      generar(A,N*N,upper);
      generar(B,M*M,upper);
/*#ifdef DEBUG
    escribir(A,N*N);
    escribir(B,M*M);
#endif*/
    ti=mseconds(); 
    printf("%d\n",sec(N,A,M,B));
    tf=mseconds(); 
      if(i!=0) tt+=tf-ti;

      free(A);
      free(B);
  }
  
    // fprintf(stats_file, "%Ld\n", tt);
    // fclose(stats_file);
  printf("%Ld\n", tt);
  return 0;
}