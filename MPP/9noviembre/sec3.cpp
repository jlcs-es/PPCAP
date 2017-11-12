/*
  CPP_CONTEST=2017
  CPP_PROBLEM=J
  CPP_LANG=XEONPHI
  CPP_PROCESSES_PER_NODE=venus 1
*/

#include <stdlib.h>
#include <omp.h>


int  sec(int n,char *a,int m,char *b)
{
    int i, j, k, h;
    int maximum=0;
    int threads;
    int size = n-m;
    int* temp =(int*) malloc(sizeof(int)*size*size); 

    #pragma offload target(mic) in(a:length(n*n)) in(b:length(m*m)) inout(temp:length(size*size)) 
    {
        threads = 4*56; // 4 x (ncore-phi - 1)
        omp_set_num_threads (threads);
        #pragma omp parallel for private(i,j,k,h)
        for (i = 0; i <= n-m; i++){
            for (j = 0; j <= n-m; j++) {
                char *aa = &a[i*n+j];
                int value=0;
                for(k=0;k < m;k++)
                    for(h=0;h < m;h++)
                        if(aa[k*n+h]==b[k*m+h])
                            value++;
                temp[i*size+j]=value;
            }
        }
    }


    //Once we have the results for each comparition we only have to know which is the best. We do this in sequencial mode.
    maximum = temp[0];
    for(int i=1; i<size*size;i++) {
    if(temp[i]>maximum)
        maximum=temp[i];
    }

    free(temp);
    return maximum;

}
