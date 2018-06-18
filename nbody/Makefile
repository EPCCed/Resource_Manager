# Compilers
CC=gcc
MCC=mcc
MPICC=mpicc
MPI_WRAPPERS=I_MPI_CC=$(MCC) MPICH_CC=$(MCC) OMPI_CC=$(MCC)

CPPFLAGS=
CFLAGS=
MCCFLAGS=
LDFLAGS=

# Nbody parameters
BIGO?=N2
BS?=2048

# Preprocessor flags
CPPFLAGS+=-Isrc -DBIGO=$(BIGO) -D_BIGO_$(BIGO) -DBLOCK_SIZE=$(BS)

# Compiler flags
CFLAGS+=-O3 -std=gnu11
MCCFLAGS+=--ompss-2 $(CFLAGS) --Wn,-O3,-std=gnu11

# Linker flags
LDFLAGS+=-lrt -lm -lstdc++

# Interoperability flags
INTEROP_DEFINE=-DINTEROPERABILITY
MPI_INTEROP_LDFLAGS=-L$(MPI_INTEROP_HOME)/lib -L$(MPI_INTEROP_HOME)/lib64 -l:libmpiompss-interop.a

# GASPI Flags
GASPI_DEVICE?=ib
GASPI_CFLAGS=-I$(GASPI_HOME)/include
GASPI_LDFLAGS=-L$(GASPI_HOME)/lib -L$(GASPI_HOME)/lib64 -lGPI2 -lpthread -lrt
GASPI_INTEROP_LDFLAGS=-L$(GASPI_INTEROP_HOME)/lib -L$(GASPI_INTEROP_HOME)/lib64 -lGPI2-ompss -lpthread -lrt

ifeq ($(GASPI_DEVICE),ib)
GASPI_LDFLAGS+=-libverbs
GASPI_INTEROP_LDFLAGS+=-libverbs
endif

PLAIN_SOURCES=          \
    src/common/common.c \
    src/plain/utils.c   \
    src/plain/main.c

SMP_SOURCES=                           \
    src/common/common.c                \
    src/blocking/common/common_utils.c \
    src/blocking/smp/utils.c           \
    src/blocking/smp/main.c

MPI_SOURCES=                           \
    src/common/common.c                \
    src/blocking/common/common_utils.c \
    src/blocking/mpi/utils.c           \
    src/blocking/mpi/main.c

GASPI_SOURCES=                         \
    src/common/common.c                \
    src/blocking/common/common_utils.c \
    src/blocking/gaspi/utils.c         \
    src/blocking/gaspi/main.c

PROGS=                                \
    nbody_seq_plain.$(BIGO).exe       \
    nbody_omp_plain.$(BIGO).exe       \
    nbody_seq.$(BIGO).$(BS).exe       \
    nbody_ompss.$(BIGO).$(BS).exe     \
    nbody_mpi.$(BIGO).$(BS).exe       \
    nbody_mpi_ompss.$(BIGO).$(BS).exe

ifdef MPI_INTEROP_HOME
PROGS+=nbody_mpi_ompss_interop.$(BIGO).$(BS).exe
endif

ifdef GASPI_HOME
PROGS+=nbody_gaspi.$(BIGO).$(BS).exe
PROGS+=nbody_gaspi_ompss.$(BIGO).$(BS).exe
endif

ifdef GASPI_INTEROP_HOME
PROGS+=nbody_gaspi_ompss_interop.$(BIGO).$(BS).exe
endif

all: $(PROGS)

nbody_seq_plain.$(BIGO).exe: $(PLAIN_SOURCES) src/plain/solver_seq.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

nbody_omp_plain.$(BIGO).exe: $(PLAIN_SOURCES) src/plain/solver_omp.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -fopenmp -o $@ $^ $(LDFLAGS)

nbody_seq.$(BIGO).$(BS).exe: $(SMP_SOURCES) src/blocking/smp/solver_seq.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

nbody_ompss.$(BIGO).$(BS).exe: $(SMP_SOURCES) src/blocking/smp/solver_ompss.c
	$(MCC) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

nbody_mpi.$(BIGO).$(BS).exe: $(MPI_SOURCES) src/blocking/mpi/solver_mpi.c
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

nbody_mpi_ompss.$(BIGO).$(BS).exe: $(MPI_SOURCES) src/blocking/mpi/solver_mpi_ompss.c
	$(MPI_WRAPPERS) $(MPICC) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

nbody_mpi_ompss_interop.$(BIGO).$(BS).exe: $(MPI_SOURCES) src/blocking/mpi/solver_mpi_ompss.c
	$(MPI_WRAPPERS) $(MPICC) $(INTEROP_DEFINE) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_INTEROP_LDFLAGS)

nbody_gaspi.$(BIGO).$(BS).exe: $(GASPI_SOURCES) src/blocking/gaspi/solver_gaspi.c
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(GASPI_CFLAGS) -o $@ $^ $(LDFLAGS) $(GASPI_LDFLAGS)

nbody_gaspi_ompss.$(BIGO).$(BS).exe: $(GASPI_SOURCES) src/blocking/gaspi/solver_gaspi_ompss.c
	$(MPI_WRAPPERS) $(MPICC) $(CPPFLAGS) $(MCCFLAGS) $(GASPI_CFLAGS) -o $@ $^ $(LDFLAGS) $(GASPI_LDFLAGS)

nbody_gaspi_ompss_interop.$(BIGO).$(BS).exe: $(GASPI_SOURCES) src/blocking/gaspi/solver_gaspi_ompss.c
	$(MPI_WRAPPERS) $(MPICC) $(INTEROP_DEFINE) $(CPPFLAGS) $(MCCFLAGS) $(GASPI_CFLAGS) -o $@ $^ $(LDFLAGS) $(GASPI_INTEROP_LDFLAGS)

check: all
	@./scripts/run-tests.sh

clean:
	rm -f *.o *.exe
