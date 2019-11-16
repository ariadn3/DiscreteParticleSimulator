#!/bin/bash
mkdir test_in
mkdir test_out
mkdir machineFiles
python3 generate.py
python3 generateMachineFile.py
mpicc -o loc source/*.c
PROCESSES=(20 40 80 160 28 56 112 224)
TOTAL_RUNS=1

for process in ${PROCESSES[@]}
do
	{
		read
		while IFS= read -r line
		do
			ssh -n "$line" -q "cd ~; mkdir -p mpiTestXeS"
			scp -r -q -o LogLevel=QUIET loc "$line":~/mpiTestXeS
		done
	} < machineFiles/"$process".mf
	for inputFile in test_in/*.in
	do
		fileName="${inputFile##*/}"
		truncFilepath=${fileName%.in}
		runTrial=1
		while [ $runTrial -le $TOTAL_RUNS ]
		do
	        	( time mpirun -machinefile machineFiles/"$process".mf -np "$process" ./loc < $inputFile > /dev/null ) 2> test_out/"$truncFilepath"-"$process"-trial"$runTrial".time
			((runTrial++))
		done
	done
done
parentdir="${PWD##*/}"
scp -r test_out e0191783@sunfire:~/testBin/"$parentdir"-"$HOSTNAME" 
cd ..
rm -rf "$parentdir"
