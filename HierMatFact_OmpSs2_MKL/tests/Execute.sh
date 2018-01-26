# Tests with 1 core (sequential)
echo "#####################################################################"	> Tests_results.txt
echo "TESTS WITH 1 CORE"							>> Tests_results.txt
echo "#####################################################################"	>> Tests_results.txt
echo "./a.out H-Matrices/matrix_600_3_00_0.txt"					>> Tests_results.txt
echo "---------------------------------------------------------------------"	>> Tests_results.txt
taskset -c 0-0 ./../a.out ./../H-Matrices/matrix_600_3_00_0.txt 		>> Tests_results.txt
echo "---------------------------------------------------------------------"	>> Tests_results.txt
echo "./a.out H-Matrices/matrix_600_3_25_0.txt"					>> Tests_results.txt
echo "---------------------------------------------------------------------"    >> Tests_results.txt
taskset -c 0-0 ./../a.out ./../H-Matrices/matrix_600_3_25_0.txt 		>> Tests_results.txt

# Tests with 2 cores
echo "#####################################################################"    >> Tests_results.txt
echo "TESTS WITH 2 CORES"                                                       >> Tests_results.txt
echo "#####################################################################"    >> Tests_results.txt
echo "./a.out H-Matrices/matrix_600_3_00_0.txt"                                 >> Tests_results.txt
echo "---------------------------------------------------------------------"    >> Tests_results.txt
taskset -c 0-1 ./../a.out ./../H-Matrices/matrix_600_3_00_0.txt			>> Tests_results.txt
echo "---------------------------------------------------------------------"    >> Tests_results.txt
echo "./a.out H-Matrices/matrix_600_3_25_0.txt"                                 >> Tests_results.txt
echo "---------------------------------------------------------------------"    >> Tests_results.txt
taskset -c 0-1 ./../a.out ./../H-Matrices/matrix_600_3_25_0.txt			>> Tests_results.txt

# Tests with 4 cores
echo "#####################################################################"    >> Tests_results.txt
echo "TESTS WITH 4 CORES"                                                       >> Tests_results.txt
echo "#####################################################################"    >> Tests_results.txt
echo "./a.out H-Matrices/matrix_600_3_00_0.txt"                                 >> Tests_results.txt
echo "---------------------------------------------------------------------"    >> Tests_results.txt
taskset -c 0-3 ./../a.out ./../H-Matrices/matrix_600_3_00_0.txt			>> Tests_results.txt
echo "---------------------------------------------------------------------"    >> Tests_results.txt
echo "./a.out H-Matrices/matrix_600_3_25_0.txt"                                 >> Tests_results.txt
echo "---------------------------------------------------------------------"    >> Tests_results.txt
taskset -c 0-3 ./../a.out ./../H-Matrices/matrix_600_3_25_0.txt			>> Tests_results.txt

echo "#####################################################################"    >> Tests_results.txt
echo "END OF TESTS"                                                       	>> Tests_results.txt
echo "#####################################################################"    >> Tests_results.txt

# Delete auxiliar files created to verify results
rm A_values.txt LU_values.txt

