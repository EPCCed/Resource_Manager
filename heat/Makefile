# Compilers
CXX=g++
MCXX=mcxx
MPICXX=mpicxx

# Use the Intel MPI compiler if possible
ifdef I_MPI_ROOT
MPICXX=mpiicpc
endif

# Set the default block size (elements)
BSX?=1024
BSY?=$(BSX)

# Preprocessor flags
CPPFLAGS=-Isrc -DBSX=$(BSX) -DBSY=$(BSY)

# Compiler flags
CFLAGS=-O3 -std=c++11
MCCFLAGS=--ompss-2 $(CFLAGS) --Wn,-O3,-std=c++11

# Linker flags
LDFLAGS=-lrt -lm

# Link to the multithreaded MPI
MPI_LDFLAGS=
ifdef I_MPI_ROOT
MPI_LDFLAGS+=-link_mpi=opt_mt
endif

# MPI Wrappers
WRAPPERS=I_MPI_CXX=$(MCXX) MPICH_CXX=$(MCXX) OMPI_CXX=$(MCXX)

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

# Extension name
ifeq ($(BSX),$(BSY))
EXT=$(BSX)bs.exe
else
EXT=$(BSX)x$(BSY)bs.exe
endif

# List of programs
PROGS=heat_seq.$(EXT)    \
    heat_ompss.$(EXT)    \
    heat_mpi.pure.$(EXT) \
    heat_mpi.omp.$(EXT)  \
    heat_mpi.task.$(EXT) \

ifdef MPI_INTEROP_HOME
PROGS+=heat_mpi.interop.$(EXT)
endif

ifdef GASPI_HOME
PROGS+=heat_gaspi.pure.$(EXT)
PROGS+=heat_gaspi.omp.$(EXT)
PROGS+=heat_gaspi.task.$(EXT)
endif

ifdef GASPI_INTEROP_HOME
PROGS+=heat_gaspi.interop.$(EXT)
endif

# Sources
SMP_SRC=src/common/misc.cpp src/smp/main.cpp
MPI_SRC=src/common/misc.cpp src/mpi/main.cpp
GASPI_SRC=src/common/misc.cpp src/gaspi/utils.cpp src/gaspi/main.cpp

all: $(PROGS)

heat_seq.$(EXT): $(SMP_SRC) src/smp/solver_seq.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

heat_ompss.$(EXT): $(SMP_SRC) src/smp/solver_ompss.cpp
	$(MCXX) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

heat_mpi.pure.$(EXT): $(MPI_SRC) src/mpi/solver_pure.cpp
	$(MPICXX) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

heat_mpi.omp.$(EXT): $(MPI_SRC) src/mpi/solver_omp.cpp
	$(WRAPPERS) $(MPICXX) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

heat_mpi.task.$(EXT): $(MPI_SRC) src/mpi/solver_task.cpp
	$(WRAPPERS) $(MPICXX) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

heat_mpi.interop.$(EXT): $(MPI_SRC) src/mpi/solver_task.cpp
	$(WRAPPERS) $(MPICXX) $(INTEROP_DEFINE) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS) $(MPI_INTEROP_LDFLAGS)

heat_gaspi.pure.$(EXT): $(GASPI_SRC) src/gaspi/solver_pure.cpp
	$(MPICXX) $(CPPFLAGS) $(CFLAGS) $(GASPI_CFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS) $(GASPI_LDFLAGS)

heat_gaspi.omp.$(EXT): $(GASPI_SRC) src/gaspi/solver_omp.cpp
	$(WRAPPERS) $(MPICXX) $(CPPFLAGS) $(MCCFLAGS) $(GASPI_CFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS) $(GASPI_LDFLAGS)

heat_gaspi.task.$(EXT): $(GASPI_SRC) src/gaspi/solver_task.cpp
	$(WRAPPERS) $(MPICXX) $(CPPFLAGS) $(MCCFLAGS) $(GASPI_CFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS) $(GASPI_LDFLAGS)

heat_gaspi.interop.$(EXT): $(GASPI_SRC) src/gaspi/solver_task.cpp
	$(WRAPPERS) $(MPICXX) $(INTEROP_DEFINE) $(CPPFLAGS) $(MCCFLAGS) $(GASPI_CFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS) $(GASPI_INTEROP_LDFLAGS)

check: all
	@./scripts/run-tests.sh

clean:
	rm -f *.o *.exe
