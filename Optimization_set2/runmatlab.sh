#!/bin/bash

# Number of files to process
num_files=20

rm -fR run-* # Careful to this line, it overwrites existing folders with this name

# Loop over files
for i in $(seq 1 $num_files); do
    mkdir "run-${i}"
    # cp cell_full.m "run-${i}"
    # cp chebdif.m "run-${i}"
    # cp cleft_bc_coupling.m "run-${i}"
    # cp cleft_odes_coupling.m "run-${i}"
    # cp clenshaw_curtis_p.m "run-${i}"
    # cp in_out.m "run-${i}"
    # cp pars5.m "run-${i}"
    # cp main_in_out.m "run-${i}"
    # cp solve_end_cleft.m "run-${i}"
    # cp objectiveFunction.m "run-${i}"

    cp *.m "run-${i}"

    cd  "run-${i}"
    # Run MATLAB code in the background
    matlab -nodisplay -nosplash -r "main_in_out; exit" > "matlab_output.log" 2>&1 &
    cd ..
done

# Wait for all background jobs to finish
wait
