#!/bin/bash
mkdir test_in
mkdir test_out
python3 generate.py
gcc source/*.h
gcc source/*.c -lm
TOTAL_RUNS=5

allThreads=(1 2 4 6 7 8 9 10 12 16 19 20 21 24 28 32 40 64)

for numThreads in ${allThreads[@]}
do
	for inputFile in test_in/*.in
	do
		fileName="${inputFile##*/}"
		truncFilepath=${fileName%.in}
		# echo "Running GPU $truncFilepath..."
		runTrial=1
		while [ $runTrial -le $TOTAL_RUNS ]
		do
			# echo -e "\tTrial $runTrial..."
			perf stat -o test_out/"$truncFilepath-$numThreads-trial$runTrial" ./a.out $numThreads < $inputFile > /dev/null
			((runTrial++))
		done
	done
done
parentdir="${PWD##*/}"
scp -r test_out e0191783@sunfire:~/testBin/"$parentdir"-"$HOSTNAME" 
cd ..
rm -rf "$parentdir"
