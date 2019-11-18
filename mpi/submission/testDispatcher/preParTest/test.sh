TOTAL_RUNS=3

allThreads=(20)

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
