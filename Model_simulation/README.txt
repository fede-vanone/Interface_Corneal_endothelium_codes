In this folder there are the codes to run the model simulation. 

The main program is 'main_coupling.m'. It produces the solution for the coupled model
of ion and water transport in the cell and in the cleft. It uses an iterative 
procedure alternating the solution for the cellular domain ('cell_full.m') 
and the solution for the cleft ('solve_end_cleft.m') until convergence. 
It produces the values for the model variables, ion fluxes and water fluxes 
and it diplays some plots.

The parameters are set in 'pars5.m'.

'solve_end_cleft.m' uses bvp5c to compute the solution in the cleft domain, 
provided the values for the cell variavles, the expression for the cleft
ODEs (in 'cleft_odes_coupling.m') and for the boundary conditions (in 
'cleft_bc_coupling.m').

'chebdif.m' and 'clenshaw_curtis_p.m' are used to approzimate integrals. 