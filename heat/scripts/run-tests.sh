#!/bin/bash

function usage {
	if [ "$#" -eq 1 ]; then
		echo "Error: $1"
	fi
	echo "Usage: ./run-tests.sh [-p exec_path] [-np num_hybrid_procs] [-nt total_threads] [-s matrix_size] [-t timesteps] [-h]"
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
size=$((8*1024))
timesteps=100

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
echo "  $size x $size matrix"
echo "  $timesteps timesteps"
echo "Getting all '*.exe' files from path '$path'..."

# Get all executable files from the desired path
programs=("$path"/*.exe)

export NANOS6=optimized

reffile=.heat_ref.ppm
tmpfile=.heat_tmp.ppm

rm -f $reffile
rm -f $tmpfile

total=0
failed=0

seq_prog=false
for ((i = 0; i < ${#programs[@]}; i++)); do
	prog="${programs[$i]}"
	progname=$(getprogname "$prog")
	if [[ $progname == *"_seq"* ]] && [ -f "$prog" ]; then
		echo "Generating reference file from ${progname}..."
		"$prog" -s $size -t $timesteps -o$reffile > /dev/null
		error=$?
		if [ "$error" -ne 0 ]; then
			echo "Error: Sequential program has failed!"
			exit 1
		fi
		seq_prog=true
	fi
done

if [ "$seq_prog" == false ]; then
	echo "Error: Sequential program not found!"
	exit 1
fi

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
	if [[ $progname == *"_mpi.pure"* ]] || [[ $progname == *"_gaspi.pure"* ]]; then
		mpiexec.hydra -n $nthreads "$prog" -s $size -t $timesteps -o$tmpfile > /dev/null
		error=$?
	elif [[ $progname == *"_mpi"* ]] || [[ $progname == *"_gaspi"* ]]; then
		mpiexec.hydra -n $nprocs -bind-to hwthread:$nthreadsxproc "$prog" -s $size -t $timesteps -o$tmpfile > /dev/null
		error=$?
	else
		"$prog" -s $size -t $timesteps -o$tmpfile > /dev/null
		error=$?
	fi

	# Check the return value of the program
	if [ "$error" -ne 0 ]; then
		failed=$(($failed + 1))
		echo $progname FAILED
		rm -f $tmpfile
		continue;
	fi

	# Result checking
	diff $reffile $tmpfile > /dev/null
	error=$?
	if [ "$error" -eq 0 ]; then
		echo $progname PASSED
	else
		failed=$(($failed + 1))
		echo $progname FAILED
	fi
	rm -f $tmpfile
done

rm -f $reffile
rm -f $tmpfile

echo ---------------------------------
if [ $failed -eq 0 ]; then
	echo SUMMARY: ALL TESTS PASSED
	exit 0
else
	echo SUMMARY: ${failed}/${total} TESTS FAILED
	exit 1
fi
echo ---------------------------------

