# N-body

## Introduction
An N-body simulation numerically approximates the evolution of a system of
bodies in which each body continuously interacts with every other body.  A
familiar example is an astrophysical simulation in which each body represents a
galaxy or an individual star, and the bodies attract each other through the
gravitational force.

N-body simulation arises in many other computational science problems as well.
For example, protein folding is studied using N-body simulation to calculate
electrostatic and van der Waals forces. Turbulent fluid flow simulation and
global illumination computation in computer graphics are other examples of
problems that use N-body simulation.

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

The nbody application has several versions which are compiled in different 
binaries, by executing the `make` command. These versions can be blocking, 
when the particle space is divided into smaller blocks, or non-blocking, when 
it is not. They are:

  * **nbody_seq_plain**: Sequential version (non-blocking).
  * **nbody_omp_plain**: Parallel version using fork-join parallelism (non-blocking).
  * **nbody_seq**: Sequential version (blocking).
  * **nbody_ompss**: Parallel version using OmpSs tasks (blocking).
  * **nbody_mpi**: Parallel version using MPI (blocking).
  * **nbody_mpi_ompss**: Parallel version using OmpSs tasks + MPI (blocking). Communication tasks are serialized.
  * **nbody_mpi_ompss_interop**: Parallel version using OmpSs tasks + MPI + Interoperability library (blocking). *See building instructions, step 1*.
  * **nbody_gaspi**: Parallel version using GASPI (blocking).
  * **nbody_gaspi_ompss**: Parallel version using OmpSs tasks + GASPI (blocking). Communication tasks are serialized.
  * **nbody_gaspi_ompss_interop**: Parallel version using OmpSs tasks + GASPI + Interoperability library (blocking). *See building instructions, step 1*.


  The simplest way to compile this package is:

  1. Stay in N-Body root directory to recursively build all the versions.
     The N-Body MPI + OmpSs tasks + Interoperability library version is
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
     Optionally, you can use a different block size when building
     the benchmark (by default 2048). Type `make BS=MY_BLOCK_SIZE`
     in order to change this value.

  3. In addition, you can type 'make check' to check the correctness
     of the built versions. By default, MPI versions run with 2 MPI
     processes and 2 hardware threads for each process, but you can
     change these parameters when executing `scripts/run-tests.sh`
     (see the available options passing `-h`).

## Execution instructions

The binaries accept several options. The most relevant options are the number 
of total particles with `-p` (default: 16384), and the number of timesteps with 
`-t` (default: 10). More options can be seen passing the `-h` option. An example 
of execution could be:

```
$ mpiexec -n 4 -bind-to hwthread:16 nbody_mpi_ompss.N2.1024.exe -t 100 -p 8192
```

in which the application will perform 100 timesteps in 4 MPI processes with 16 
hardware threads in each process (used by the OmpSs runtime). The total number 
of particles will be 8192, this means that each process will have 2048 
particles (2 blocks per process).

