This folder contains the codes for computing the first and total order sensitivity indices.
The model we consider is the coupled model cell-cleft. 
In the subfolder Analysis there are the codes for analysing the results and displaying plots.

1) Set the parameters

Parameter_settings_EFAST:
Sets all the parameters in the model, including those you want to span. 
pmean - is the mean value of the parameters
pmin - minimum value
pmax - maximum value

selectX - spans parameters with given frequencies and saves them into the variable X and the file X.mat.
It uses functions SETFREQ.m to set the frequencies and parameterdist.m to set 
distribution of the parameters (you can pick between uniform, normal or lognormal).

2) Run the solver for parameters from X.mat

This is the file XtoY.m, which calls the function model_main.m. 
The function model_main.m runs the cellular model coupled_model.m with the permeabilities 
from matrix X corresponding to the run_num specified as input. The other parameters 
needed for model_main.m are taken from pars5.m. cell_full.m contains the cellular model,
cleft_bc_coupling.m the boundary conditions for the cleft ODEs, which are 
specified in cleft_odes_coupling.m. 
The output is stored in variable Y and the file Yfull.mat. You can pick multiple outputs,
like various concentrations, fluxes, etc.
To make things quicker you can run in parallel for each parameter (variable i=1:k).
To run on the background you can run the script run_XtoY from the terminal.
chebdif.m and clenshaw_curtis_p2.m are used for integration. 

3) Analyse the results

Once you have your output Y, you can analyse it. Here Yfull.mat is the output from XtoY.m
and X.mat is the parameter file. 
To run SA, you can use figure_wat_graf.m file. It produces bar plots with sensitivity indices. 
These indices are calculated by function efast_sd_2.m 
interp_nonan.m takes care of NANs if you have any: interp_nonan.m.
Finally, to get scatter plots for each  parameter you can run script some_plots.m.
The last script also displays scatterplots with polynomial fitting for TEP vs
water flux, curent across the tight junction vs water flux and current across the
tight junction vs TEP.
