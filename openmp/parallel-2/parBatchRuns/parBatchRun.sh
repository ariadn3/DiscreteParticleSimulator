#!/bin/bash

python3 generateTestcases.py

TOTAL_RUNS=5

rm *.par -f

allThreads=(1 2 4 6 7 8 9 10 12 16 19 20 21 24 28 32 40 64)
gridSizes=(1 2 4 8 16 32 64)

for numThreads in ${allThreads[@]}
do
	echo "Running on $numThreads threads..."
	for gSize in ${gridSizes[@]}
	do
		echo -e "\tRunning on grid size $gSize..."
		for inputFile in *.in
		do
			truncFilepath=${inputFile%.in}
			echo -e "\t\tRunning parallel $truncFilepath..."
			runTrial=1
			while [ $runTrial -le $TOTAL_RUNS ]
			do
				echo -e "\t\t\tTrial $runTrial..."
				perf stat -o "$truncFilepath-$numThreads-$gSize-trial$runTrial.par" ../a.out $numThreads $gSize < $inputFile > /dev/null
				((runTrial++))
			done
		done
	done
done