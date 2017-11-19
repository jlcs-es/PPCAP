CC = icc
C_FLAGS	= -O3 -qopenmp

EXE = matriz_matriz_tras
LD_LIBS	= -qoffload-option,mic,ld, "-L$(MKLROOT)/lib/mic -Wl, \
		  --start-group -lmkl_intel_lp64 -lmkl_intel_thread \
		  -lmkl_core -Wl,--end-group" -mkl -lm

$(EXE): *.o 
	$(CC) $^ -o $@ $(LD_LIBS)

%.o: %.c
	$(CC) $(C_FLAGS) -c $^ -o $@

clean : 
	rm -f *.o $(EXE)