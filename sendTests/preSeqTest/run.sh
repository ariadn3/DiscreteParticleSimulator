#!/bin/bash
mkdir test_in
mkdir test_out
python3 generate.py
gcc source/*.h
gcc source/*.c -lm
TOTAL_RUNS=5

for inputFile in test_in/*.in
do
	fileName="${inputFile##*/}"
	truncFilepath=${fileName%.in}
	# echo "Running GPU $truncFilepath..."
	runTrial=1
	while [ $runTrial -le $TOTAL_RUNS ]
	do
		# echo -e "\tTrial $runTrial..."
		perf stat -o test_out/"$truncFilepath-trial$runTrial" ./a.out < $inputFile > /dev/null
		((runTrial++))
	done
done
parentdir="${PWD##*/}"
scp -r test_out e0191783@sunfire:~/testBin/"$parentdir"-"$HOSTNAME" 
cd ..
rm -rf "$parentdir"
