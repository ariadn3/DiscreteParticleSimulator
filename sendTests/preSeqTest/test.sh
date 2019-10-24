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
