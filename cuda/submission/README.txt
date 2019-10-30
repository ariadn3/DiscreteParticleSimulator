To execute the code:

1.  Proceed to either /code/cpu or /code/gu

2.  (For OpenMP implementation) Compile the C code with the following commands

    (i)  gcc *.h *.c -fopenmp -lm

3.  (For CUDA implementation) Compile the C code with the following commands

    (i)  (FMA enabled) nvcc *.cu -o fast.exe -arch sm_75 --compiler-options -Wall --fmad=true
    (ii) (FMA DISABLED) nvcc *.cu -o fast.exe -arch sm_75 --compiler-options -Wall --fmad=false