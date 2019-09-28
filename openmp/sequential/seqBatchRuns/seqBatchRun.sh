#!/bin/bash

TOTAL_RUNS=5

rm *.seq -f
for inputFile in *.in
do
	truncFilepath=${inputFile%.in}
	echo "Running sequential $truncFilepath..."
	runTrial=1
	while [ $runTrial -le $TOTAL_RUNS ]
	do
		echo -e "\tTrial $runTrial..."
		perf stat -o "$truncFilepath-trial$runTrial.seq" ../a.out < $inputFile > /dev/null
		((runTrial++))
	done
done
