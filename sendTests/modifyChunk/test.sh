TOTAL_RUNS=5

oneDChunkSizes=(32 64 96 128 192 256)
twoDChunkSizes=(32 64 96 128 192 256)

for oneDChunkSize in ${oneDChunkSizes[@]}
do
    for twoDChunkSize in ${twoDChunkSizes[@]}
    do
        for inputFile in test_in/*.in
        do
	        fileName="${inputFile##*/}"
	        truncFilepath=${fileName%.in}
	        runTrial=1
	        while [ $runTrial -le $TOTAL_RUNS ]
	        do
	        	/usr/local/cuda/bin/nvprof --csv --log-file test_out/"$truncFilepath"-"$oneDChunkSize"-"$twoDChunkSize"-fmad-trial"$runTrial" ./withfMad < $inputFile > /dev/null
        		/usr/local/cuda/bin/nvprof --csv --log-file test_out/"$truncFilepath"-"$oneDChunkSize"-"$twoDChunkSIze"-nofmad-trial"$runTrial" ./withoutfMad < $inputFile > /dev/null
     		((runTrial++))
     	    done
        done
    done
done
