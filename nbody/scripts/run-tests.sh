#!/bin/bash

function usage {
	if [ "$#" -eq 1 ]; then
		echo "Error: $1"
	fi
	echo "Usage: ./run-tests.sh [-p exec_path] [-np num_hybrid_procs] [-nt total_threads] [-s num_particles] [-t timesteps] [-h]"
	exit 1
}

function getprogname {
	local progarray=()
	IFS='/' read -r -a progarray <<< "$1"
	local index=$((${#progarray[@]} - 1))
	echo "${progarray[$index]}"
}

# Default executables path
path=$PWD

# Default configuration
nprocs=2
nthreads=4

# Default execution parameters
size=$((16*1024))
timesteps=20

nargs=$#
args=("$@")

error=$(($nargs % 2))
if [ "$error" -eq 1 ]; then
	usage
fi

i=0
while (("$i" < "$nargs")); do
	opt=${args[$i]}
	((i++))
	val=${args[$i]}
	((i++))

	if [ "$opt" == "-p" ]; then
		path=$val
	elif [ "$opt" == "-np" ]; then
		nprocs=$val
	elif [ "$opt" == "-nt" ]; then
		nthreads=$val
	elif [ "$opt" == "-s" ]; then
		size=$val
	elif [ "$opt" == "-t" ]; then
		timesteps=$val
	else
		usage "'$opt' is not a valid option!"
	fi
done

nthreadsxproc=$((nthreads / nprocs))

echo ---------------------------------
echo TEST SUITE
echo ---------------------------------

echo "Configuration (Distributed versions):"
echo "  $nthreads processes"
echo "Configuration (Hybrid versions):"
echo "  $nprocs processes"
echo "  $nthreadsxproc threads per proc"
echo "Execution parameters:"
echo "  $size particles"
echo "  $timesteps timesteps"
echo "Getting all '*.exe' files from path '$path'..."

# Get all executable files from the desired path
programs=("$path"/*.exe)

# Oss loop is only supported by "fifo" and "naive" schedulers
export NANOS6_SCHEDULER=fifo
export NANOS6=optimized

logfile=.nbody_log

rm -f $logfile

total=0
failed=0

echo "Starting the verification..."
echo ---------------------------------

for ((i = 0; i < ${#programs[@]}; i++)); do
	prog="${programs[$i]}"
	# Check if it exists
	if [ ! -f "$prog" ]; then
		continue;
	fi
	total=$(($total + 1))

	progname=$(getprogname "$prog")

	# Execution
	if [[ $progname == *"_mpi."* ]] || [[ $progname == *"_gaspi."* ]]; then
		mpiexec.hydra -n $nthreads "$prog" -p $size -t $timesteps -c &> $logfile
		error=$?
	elif [[ $progname == *"_mpi_ompss"* ]] || [[ $progname == *"_gaspi_ompss"* ]]; then
		mpiexec.hydra -n $nprocs -bind-to hwthread:$nthreadsxproc "$prog" -p $size -t $timesteps -c &> $logfile
		error=$?
	else
		"$prog" -p $size -t $timesteps -c &> $logfile
		error=$?
	fi

	# Check the return value of the program
	if [ "$error" -ne 0 ]; then
		failed=$(($failed + 1))
		echo $progname FAILED
		rm -f $logfile
		continue;
	fi

	# Result checking
	grep -Fxq 'Result validation: OK' $logfile
	error=$?
	if [ "$error" -eq 0 ]; then
		echo $progname PASSED
	else
		failed=$(($failed + 1))
		echo $progname FAILED
	fi
	rm -f $logfile
done

rm -f $logfile

echo ---------------------------------
if [ $failed -eq 0 ]; then
	echo SUMMARY: ALL TESTS PASSED
	exit 0
else
	echo SUMMARY: ${failed}/${total} TESTS FAILED
	exit 1
fi
echo ---------------------------------
