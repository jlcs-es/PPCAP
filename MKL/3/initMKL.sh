#!/bin/bash

source /etc/profile.d/modules.sh
module load intel/2017_u1
source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64
#module load intel

# icc -O3 code.c -mkl

# icc -O3 lu_mkl_offload.c io.c -I/opt/intel/composer_xe_2015.3.187/mkl/include -qoffloadoption,mic,ld,"-L$(MKLROOT)/lib/mic -Wl, --start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group" -mkl
icc -O3 lu_mkl_offload.c io.c -qoffload-option,mic,ld, "-L$MKLROOT/lib/mic -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group" -mkl -lm

# export LD_LIBRARY_PATH
# export MIC_LD_LIBRARY_PATH
# export MIC_ENV_PREFIX=PHI
# export PHI_KMP_AFFINITY=granularity=fine,balanced
