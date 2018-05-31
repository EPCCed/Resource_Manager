STARPU_HOME = $(HOME)/binaries/nord/starpu
NANOS6_HOME = $(HOME)/binaries/interop-nord/nanos6
OPENCL_HOME = $(HOME)/binaries/nord/opencl
HWLOC_HOME = $(HOME)/binaries/nord/hwloc

LASS_HOME = LASs/
LASS_FLAGS = -DLASs_WITH_MKL

OUT_DIR = ./bin
SRC_DIR = ./src
FLAGS = -O3 -Wall $(LASS_FLAGS)

.PHONY: all outdir lass clean check

all: outdir lass $(OUT_DIR)/interop_starpu_ompss

opencl: outdir lass $(OUT_DIR)/interop_starpu_opencl_ompss

$(OUT_DIR)/interop_starpu_opencl_ompss: $(OUT_DIR)/ompss_dgemm.o $(OUT_DIR)/opencl_dgemm.o $(OUT_DIR)/starpu_opencl_dgemm.o $(OUT_DIR)/main.o
	$(CC) $(FLAGS) -L$(STARPU_HOME)/lib -L$(MKL_HOME)/lib/intel64 -L$(LASS_HOME) -L$(NANOS6_HOME)/lib -L$(OPENCL_HOME)/lib -L$(HWLOC_HOME)/lib $? -o $@ -Wl,-Bstatic -llass -Wl,-Bdynamic -lstarpu-1.3 -lstarpurm-1.3 -lmkl_sequential -lmkl_core -lmkl_rt -lnanos6 -lpocl

$(OUT_DIR)/interop_starpu_ompss: $(OUT_DIR)/ompss_dgemm.o $(OUT_DIR)/starpu_dgemm.o $(OUT_DIR)/main.o
	$(CC) $(FLAGS) -DENABLE_OPENCL -L$(STARPU_HOME)/lib -L$(MKL_HOME)/lib/intel64 -L$(LASS_HOME) -L$(NANOS6_HOME)/lib -L$(HWLOC_HOME)/lib $? -o $@ -Wl,-Bstatic -llass -Wl,-Bdynamic -lstarpu-1.3 -lstarpurm-1.3 -lmkl_sequential -lmkl_core -lmkl_rt -lnanos6

$(OUT_DIR)/ompss_dgemm.o: $(SRC_DIR)/ompss_dgemm.c
	$(CC) $(FLAGS) -I$(LASS_HOME)/include -I$(NANOS6_HOME)/include -I$(MKL_HOME)/include $? -c -o $@

$(OUT_DIR)/starpu_opencl_dgemm.o: $(SRC_DIR)/starpu_dgemm.c
	$(CC) $(FLAGS) -I$(STARPU_HOME)/include/starpu/1.3 -I$(STARPU_HOME)/include/starpurm/1.3 -I$(HWLOC_HOME)/include -I$(MKL_HOME)/include $? -c -o $@

$(OUT_DIR)/starpu_dgemm.o: $(SRC_DIR)/starpu_dgemm.c
	$(CC) $(FLAGS) -I$(STARPU_HOME)/include/starpu/1.3 -I$(STARPU_HOME)/include/starpurm/1.3 -I$(HWLOC_HOME)/include -I$(MKL_HOME)/include $? -c -o $@

$(OUT_DIR)/opencl_dgemm.o: $(SRC_DIR)/opencl_dgemm.c
	$(CC) $(FLAGS) -I$(LASS_HOME)/include -I$(OPENCL_HOME)/include -I$(MKL_HOME)/include $? -c -o $@

$(OUT_DIR)/main.o: $(SRC_DIR)/main.c
	$(CC) $(FLAGS) -I$(MKL_HOME)/include $? -c -o $@

outdir: $(OUT_DIR)

lass: LASs
	$(MAKE) -C LASs

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

clean:
	rm -rf $(OUT_DIR)/*.o $(OUT_DIR)/interop_starpu_opencl_ompss

check: outdir $(OUT_DIR)/interop_starpu_ompss
	test/run_test_intertwine.sh
	
