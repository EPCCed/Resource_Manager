COMP=gcc

MKL_LIB_DIR = /mnt/volume100/intel/mkl/lib/intel64
MKL_INC_DIR = /mnt/volume100/intel/mkl/include
CFLAGS = -Wall
INCS = -I$(MKL_INC_DIR)
LIBS = -L$(MKL_LIB_DIR) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_rt -lpthread -lm

all: HLU_OpenMP

HLU_OpenMP: src/HLU_OpenMP.c 
	$(COMP) -fopenmp $(INCS) src/HLU_OpenMP.c src/Read_matrix_from_file.c src/TreeFunctions.c src/Verify.c $(LIBS)
clean:
	rm -f a.out
