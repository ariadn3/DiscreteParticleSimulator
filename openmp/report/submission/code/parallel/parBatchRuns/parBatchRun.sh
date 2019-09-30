#!/bin/bash

python3 generateTestcases.py

TOTAL_RUNS=5

rm *.seq -f

allThreads=(1 2 4 6 7 8 9 10 12 16 19 20 21 24 28 32 40 64)

for numThreads in ${allThreads[@]}
do
	echo "Running on $numThreads threads..."
	for inputFile in *.in
	do
		truncFilepath=${inputFile%.in}
		echo -e "\tRunning parallel $truncFilepath..."
		runTrial=1
		while [ $runTrial -le $TOTAL_RUNS ]
		do
			echo -e "\t\tTrial $runTrial..."
			perf stat -o "$truncFilepath-$numThreads-trial$runTrial.par" ../a.out $numThreads < $inputFile > /dev/null
			((runTrial++))
		done
	done
done
