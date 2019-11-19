To execute the code:

1.  Proceed to either /code

2. Compile the C code with the following commands

    (i) mpicc *.h *.c -lm

3. Run the code with the following commands

    (i) mpirun -machinefile machineFiles/<numOfProcesses>.mf -np <numOfProcesses> ./a.out

