PROCESSES=(8 16 32 64)
TOTAL_RUNS=5

for process in ${PROCESSES[@]}
do
	{
		read
		while IFS= read -r line
		do
			ssh -n "$line" -q "cd ~; mkdir -p mpiTestI7"
			scp -r -q -o LogLevel=QUIET a.out "$line":~/mpiTestI7
		done
	} < machineFiles/"$process".mf
	for inputFile in test_in/*.in
	do
		fileName="${inputFile##*/}"
		truncFilepath=${fileName%.in}
		runTrial=1
		while [ $runTrial -le $TOTAL_RUNS ]
		do
			perf stat -o test_out/"$truncFilepath-$process-trial$runTrial" mpirun -machinefile machineFiles/"$process".mf -np "$process" ./a.out $numThreads < $inputFile > /dev/null
			((runTrial++))
		done
	done
done
