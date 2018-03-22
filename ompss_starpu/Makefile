STARPU_HOME = $(HOME)/binaries/nord/starpu/
CHAMELEON_HOME = $(HOME)/binaries/nord/chameleon/
HWLOC_HOME =  $(HOME)/binaries/nord/hwloc/

FLAGS       = -DLASs_WITH_CHAMELEON

CC = mcc
SRC_DIR = src
OUT_DIR = bin

CFLAGS = -O3 -Wall --ompss-2 $(FLAGS) -I./include -I$(MKL_HOME)/include -I$(STARPU_HOME)/include/starpu/1.3 -I$(STARPU_HOME)/include/starpurm/1.3 -I$(CHAMELEON_HOME)/include/chameleon -I$(HWLOC_HOME)/include

LDFLAGS = -L$(MKL_HOME)/lib/intel64 -L$(STARPU_HOME)/lib -L$(CHAMELEON_HOME)/lib -L$(HWLOC_HOME)/lib -lm -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lchameleon -lchameleon_starpu -lhwloc -lstarpurm-1.3 -lstarpu-1.3 -lcoreblas

.PHONY: all outdir clean check

all: outdir $(OUT_DIR)/test_dpotrf

$(OUT_DIR)/test_dpotrf:
	#TILE TEST
	$(CC) $(CFLAGS) $(LDFLAGS) -c $(SRC_DIR)/ddss_tile.c -o $(OUT_DIR)/ddss_tile.o
	$(CC) $(CFLAGS) $(LDFLAGS) -c $(SRC_DIR)/ddss_flat2tiled.c -o $(OUT_DIR)/ddss_flat2tiled.o
	$(CC) $(CFLAGS) $(LDFLAGS) -c $(SRC_DIR)/ddss_tiled2flat.c -o $(OUT_DIR)/ddss_tiled2flat.o
	
	#INTEROPERABILITY
	$(CC) $(CFLAGS) $(LDFLAGS) -c $(SRC_DIR)/starpu_interop.c -o $(OUT_DIR)/starpu_interop.o
	$(CC) $(CFLAGS) $(LDFLAGS) -c $(SRC_DIR)/interop_cpu_mask.c -o $(OUT_DIR)/interop_cpu_mask.o
	
	#DPOTRF TEST
	$(CC) $(CFLAGS) $(LDFLAGS) -c $(SRC_DIR)/ddss_dpotrf.c -o $(OUT_DIR)/ddss_dpotrf.o
	$(CC) $(CFLAGS) $(LDFLAGS) -c $(SRC_DIR)/kdpotrf.c -o $(OUT_DIR)/kdpotrf.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(OUT_DIR)/ddss_tile.o \
		  $(OUT_DIR)/ddss_flat2tiled.o $(OUT_DIR)/ddss_tiled2flat.o \
		  $(OUT_DIR)/ddss_dpotrf.o $(OUT_DIR)/kdpotrf.o \
		  $(OUT_DIR)/starpu_interop.o $(OUT_DIR)/interop_cpu_mask.o \
		  $(SRC_DIR)/main.c -o $(OUT_DIR)/test_dpotrf

outdir: $(OUT_DIR)

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

clean:
	rm $(OUT_DIR)/*

check: outdir $(OUT_DIR)/test_dpotrf
	test/run_test_intertwine.sh
