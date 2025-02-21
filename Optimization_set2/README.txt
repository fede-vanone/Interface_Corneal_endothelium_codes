This folder contains the codes for the optimization of the full model. The set
of permeabilities obtained with the codes as in the folder is set 2.
Execute the file runmatlab.sh from the terminal. This will create 20 copies of
the folder Optimization_set2, called run-i, with i from 1 to 20. Then, it will
start one Matlab session in each folder, running 'main_in_out.m'.
'main_in_out.m' performs a constrained optimization using fmincon. As
starting point it takes a random perturbation of the vector input_ref (which is the 
vector of permeabilities of set 1).
Once optimization has finished, the optimum is stored in 'optimum_perm.mat' in
the corresponding folder, and the value of the objective function at the last 
iteration is stored in 'objective.mat'. These last values can be plotted against 
the number of the folder using the script 'check_obj_value.m'. 
To monitor the advance of the optimization in nodisplay mode, we write a 
file called 'matlab_output.log' per each folder.
We upload here the run-i folders (containing only output files) we obtained
with our codes. By running 'runmatlab.sh' the folders will be overwritten.
Set 2 correspomnds to the permeabilities found in folder run-1.
