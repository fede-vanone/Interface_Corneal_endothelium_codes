This folder contains the codes for optimizing the cellular permeabilities 
Ppump,Pnkcc,Pkirb,Pkira,Pcftr,Pn2b,Pn3b,Paea,Pae,P0tj,P1tj,P2tj,P3tj in the 
cellular domain. The set obtained with the codes as in this folder is set 1.

'main_opt_nd.m' is the main code. Running it performs the optimization using
the code 'par_opt_nd.m', which implements the optimization procedure fminsearch
with the objective function specified by 'objective.m'. 
Once the optimal set is found, the cellular equations are run with the 
optimized values to compute the value of the cellular variables and the ion
fluxes. 

'pars5.m' defines the parameters. 