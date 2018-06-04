# Heat

## Introduction
The Heat simulation uses an iterative Gauss-Seidel method to solve the heat equation,
which is a parabolic partial differential equation that describes the distribution of
heat (or variation in temperature) in a given region over time.

The heat equation is of fundamental importance in a wide range of science fields. In
mathematics, it is the parabolic partial differential equation par excellence. In statistics,
it is related to the study of the Brownian motion. Also, the diffusion equation is a generic
version of the heat equation, and it is related to the study of chemical diffusion processes.

## Requirements
The requirements of this application are shown in the following lists. The main requirements are:

  * **GNU Compiler Collection**.
  * **OmpSs-2**: OmpSs-2 is the second generation of the OmpSs programming model. It is a task-based
    programming model originated from the ideas of the OpenMP and StarSs programming models. The
    specification and user-guide are available at https://pm.bsc.es/ompss-2-docs/spec/ and
    https://pm.bsc.es/ompss-2-docs/user-guide/, respectively. OmpSs-2 requires both Mercurium and
    Nanos6 tools. Mercurium is a source-to-source compiler which provides the necessary support for
    transforming the high-level directives into a parallelized version of the application. The Nanos6
    runtime system library provides the services to manage all the parallelism in the application (e.g. task
    creation, synchronization, scheduling, etc). Both can be downloaded from https://github.com/bsc-pm.
  * **MPI**: This application requires an MPI library supporting the multi-threaded mode. It mainly targets
    MPICH and Intel MPI implementations. However, it should work with other libraries by adding the needed
    implementation-specific flags in the Makefile.

In addition, there are some optional tools which enable the building of other application versions:

  * **Task-Aware MPI (TAMPI)**: The Task-Aware MPI library provides the interoperability mechanism for MPI
    and OmpSs. Please contact <pm-tools@bsc.es> to get access to the TAMPI library.
  * **GASPI**: A GASPI implementation can also be used in this application. It mainly targets the GPI-2 library.
    However, the implementation must provide support for Queue Groups. In addition, there is an
    extended GASPI version which allows the interoperability between GASPI and OmpSs. Please contact
    <pm-tools@bsc.es> to get access to both extended implementations.

## Available versions and building instructions

The heat application has several versions which are compiled in different 
binaries, by executing the `make` command. They are:

  * **heat_seq**: Sequential version.
  * **heat_ompss**: Parallel version using OmpSs tasks.
  * **heat_mpi.pure**: Parallel version using MPI.
  * **heat_mpi.omp**: Parallel version using MPI + OmpSs. Fork-join parallelization.
  * **heat_mpi.task**: Parallel version using MPI + OmpSs tasks. Communication tasks are serialized.
  * **heat_mpi.interop**: Parallel version using MPI + OmpSs tasks + Interoperability library. *See building instructions, step 1*.
  * **heat_gaspi.pure**: Parallel version using GASPI.
  * **heat_gaspi.omp**: Parallel version using GASPI + OmpSs. Fork-Join parallelization.
  * **heat_gaspi.task**: Parallel version using GASPI + OmpSs tasks. Communication tasks are serialized.
  * **heat_gaspi.interop**: Parallel version using MPI + OmpSs tasks + Interoperability. *See building instructions, step 1*.


  The simplest way to compile this package is:

  1. Stay in Heat root directory to recursively build all the versions.
     The Heat MPI + OmpSs tasks + Interoperability library version is
     compiled only if the environment variable `MPI_INTEROP_HOME`
     is set to the Task-Aware MPI (TAMPI) library's installation directory.

     GASPI versions are compiled only if the environment variable `GASPI_HOME`
     is set to the GASPI's installation directory. In addition, the
     GASPI+OmpSs+Interoperability version is compiled only if the environment variable
     `GASPI_INTEROP_HOME` is set to the directory where resides the installation
     of GASPI version that supports the interoperability mechanism.

     By default the application is compiled for Infiniband (`GASPI_DEVICE=ib`) linking
     to the `libverbs` library. This can be changed with `make GASPI_DEVICE=tcp`.

  2. Type `make` to compile the selected benchmark's version(s).
     Optionally, you can use a different block size in each dimension
     (BSX and BSY for vertical and horizontal dimensions, respectively)
     when building the benchmark (by default 1024). Type
     `make BSX=MY_BLOCK_SIZE_X BSY=MY_BLOCK_SIZE_Y` in order to change
     this value. If you want the same value in each dimension, type
     `make BSX=MY_BLOCK_SIZE`.

  3. In addition, you can type 'make check' to check the correctness
     of the built versions. By default, the pure MPI version runs with
     4 processes and the hybrid versions run with 2 MPI processes and 2
     hardware threads for each process. You can change these parameters
     when executing `scripts/run-tests.sh` (see the available options
     passing `-h`).

## Execution instructions

The binaries accept several options. The most relevant options are the size 
of the matrix with `-s` (default: 2048), and the number of timesteps with 
`-t` (default: 100). More options can be seen passing the `-h` option. An example 
of execution could be:

```
$ mpiexec -n 4 -bind-to hwthread:16 heat_mpi.task.1024bs.exe -t 150 -s 8192
```

in which the application will perform 150 timesteps in 4 MPI processes with 16 
hardware threads in each process (used by the OmpSs runtime). The size of the
matrix in each dimension will be 8192 (8192^2 elements in total), this means
that each process will have 2048 * 8192 elements (16 blocks per process).

