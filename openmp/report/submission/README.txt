To execute the code:

1.  Proceed to either /code/sequential or /code/parallel

2.  Compile the C code with the following commands

    (i)  gcc *.h
    (ii) gcc -fopenmp *.c -lm

3.  (If executing own testcases)

    (i)  Sequential: ./a.out < [testfile]
    (ii) Parallel: ./a.out [number of threads] < [testfile]

4.  (For benchmarking with batch runs) Proceed to the respective batchRuns sub-folder

    (i)  Execute the command ./seqBatchRun.sh (sequential) or ./parBatchRun.sh (parallel)

To visualise the output from a testcase (requires 'print' simulation mode):

1.  Install the imagemagick Python module with pip3

2.  Run python3 generateAnimation.py with the input testcase file and output file in the same folder