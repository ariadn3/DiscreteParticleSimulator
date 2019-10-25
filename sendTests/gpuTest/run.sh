#!/bin/bash
mkdir test_in
mkdir test_out
python3 generate.py
/usr/local/cuda/bin/nvcc -o withfMad source/*.cu -fmad=true
/usr/local/cuda/bin/nvcc -o withoutfMad source/*.cu -fmad=false
TOTAL_RUNS=5

for inputFile in test_in/*.in
do
	fileName="${inputFile##*/}"
	truncFilepath=${fileName%.in}
	runTrial=1
	while [ $runTrial -le $TOTAL_RUNS ]
	do
		# time ./withfMad < $inputFile > test_out/"$truncFilepath"-fmad-trial"$runTrial".time > /dev/null
        ( time /usr/local/cuda/bin/nvprof --csv --log-file test_out/$truncFilepath-fmad-trial$runTrial ./withfMad < $inputFile > /dev/null ) 2> test_out/"$truncFilepath"-fmad-trial"$runTrial".time
		# time ./withoutfMad < $inputFile > test_out/"$truncFilepath"-nofmad-trial"$runTrial".time > /dev/null
		( time /usr/local/cuda/bin/nvprof --csv --log-file test_out/$truncFilepath-nofmad-trial$runTrial ./withoutfMad < $inputFile > /dev/null ) 2> test_out/"$truncFilepath"-nofmad-trial"$runTrial".time
		((runTrial++))
	done
done
parentdir="${PWD##*/}"
scp -r test_out e0191783@sunfire:~/testBin/"$parentdir"-"$HOSTNAME" 
cd ..
rm -rf "$parentdir"
