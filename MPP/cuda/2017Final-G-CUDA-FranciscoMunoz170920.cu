
/*
  CPP_CONTEST=2017
  CPP_PROBLEM=I
  CPP_LANG=CUDA
  CPP_PROCESSES_PER_NODE=saturno 1
*/



/* RECORD
Francisco Muñoz García
September 20, 2017
in CESGA
time 1520
speed-up 9.80
*/



#include <stdlib.h>

__device__ int count(int ld,int n,char *a,char *b) //Each CUDA thread do this work and is called from kernel so we change to __device__
{
  int i,j;
  int value=0;
  for(i=0;i < n;i++)
    for(j=0;j < n;j++)
      if(a[i*ld+j]==b[i*n+j])
        value++;
  return value;
}

/*
We create one thread for each element in matrix sizexsize. Each element compare its matrix and save the results in a matrix. For that reason
each thread has an associated element in the matrix.
*/
__global__ void mask(char* a, char* b, int* temp, int n, int m) {
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;
	int size = n-m;
	if((i<size) && (j<size)) {
		temp[i*size+j]=count(n,m,&a[i*n+j],b);
	}
}

int  sec(int n,char *a,int m,char *b)
{
  int i, j;
  int maximum=0,value;
  int size = n-m;
  int nbytes_a = sizeof(char)*n*n;
  int nbytes_b = sizeof(char)*m*m;
  int nBytes_temp = sizeof(int)*size*size;

  int* temp =(int*) malloc(sizeof(int)*size*size); 
  int* temp_d;
  char* a_d;
  char* b_d;

  int bl_dim1 = 4;
  int bl_dim2 = 8;

  dim3 block(bl_dim1,bl_dim2);

  //we need n-m threads

  int gsx = size / bl_dim1;
  if(size%bl_dim1) gsx++;
  int gsy = size / bl_dim2;
  if(size%bl_dim2) gsy++;
  dim3 grid(gsx, gsy);


  //We reserve memory for GPU
  cudaMalloc((void **) &temp_d, nBytes_temp);
  cudaMalloc((void**) &a_d, nbytes_a);
  cudaMalloc((void**) &b_d, nbytes_b);
  
  //Transfers here
  cudaMemset(temp_d, 0, nBytes_temp*sizeof(char)); //All the values should stat with zeros because each thread add values from that initial zero.

  cudaMemcpy(a_d, a, nbytes_a, cudaMemcpyHostToDevice);
  cudaMemcpy(b_d, b, nbytes_b, cudaMemcpyHostToDevice);



  //call the kernel
  mask<<<grid, block>>>(a_d, b_d, temp_d, n,m );


  //We transfer the results to RAM
  cudaMemcpy(temp, temp_d, nBytes_temp, cudaMemcpyDeviceToHost);

  cudaFree((void**)temp_d);
  cudaFree((void**)a_d);
  cudaFree((void**)b_d);

  //Once we have the results for each comparition we only have to know which is the best. We do this in sequencial mode.
  maximum = temp[0];
  for(int i=1; i<size*size;i++) {
	if(temp[i]>maximum)
		maximum=temp[i];
  }

  free(temp);
  return maximum;
}


