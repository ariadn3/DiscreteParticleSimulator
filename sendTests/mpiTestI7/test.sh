PROCESSES=(8 16 32 64 28 56 112 224)
TOTAL_RUNS=1

for process in ${PROCESSES[@]}
do
	{
		read
		while IFS= read -r line
		do
			ssh -n "$line" -q "cd ~; mkdir -p mpiTestI7"
			scp -r -q -o LogLevel=QUIET loc "$line":~/mpiTestI7
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
