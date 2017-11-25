CC = mpicc
C_FLAGS	= -O3 -qopenmp

OMPI_DIR = /usr/local/openmpi/1.6.4/intel

INC	= -I$(MKLROOT)/include \
	  -I$(OMPI_DIR)/include

EXE = matriz_matriz_mpi_omp
LD_LIBS	= -qoffload-option,mic,ld, "-L$(MKLROOT)/lib/mic -Wl, \
		  --start-group -lmkl_intel_lp64 -lmkl_intel_thread \
		  -lmkl_core -Wl,--end-group" -mkl -lm

$(EXE): $(EXE).o
	$(CC) $< -o $@ $(LD_LIBS)

%.o: %.c
	$(CC) $(C_FLAGS) $(INC) -c $< -o $@

clean : 
	rm -f *.o $(EXE)