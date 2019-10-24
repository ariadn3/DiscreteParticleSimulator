TOTAL_RUNS=5

for inputFile in test_in/*.in
do
	fileName="${inputFile##*/}"
	truncFilepath=${fileName%.in}
	runTrial=1
	while [ $runTrial -le $TOTAL_RUNS ]
	do
		/usr/local/cuda/bin/nvprof --csv --log-file test_out/"$truncFilepath"-fmad-trial"$runTrial" ./withfMad < $inputFile > /dev/null
		/usr/local/cuda/bin/nvprof --csv --log-file test_out/"$truncFilepath"-nofmad-trial"$runTrial" ./withoutfMad < $inputFile > /dev/null
		((runTrial++))
	done
done
