#!/bin/bash

echo "---------------------------------"
echo "TEST SUITE"
echo "---------------------------------"

passed=0
total=1

path=$(dirname $0)
bin_name=test_dpotrf
$path/../bin/${bin_name} 8192 U >& /dev/null

if [ $? == 0 ]; then
	echo "$bin_name PASSED"
	passed=$(($passed + 1))
else
	echo "$bin_name FAILED"
fi

echo "---------------------------------"
if [ $passed != $total ]; then
	echo "SUMMARY: ${passed}/${total} TESTS FAILED"
else
	echo "SUMMARY: ALL TESTS PASSED"
fi
echo "---------------------------------"

