function par_vec=par_opt_nd(par_guess)
% this function runs fminsearch for the function 'objective.m' and displays
% the plot of the optimization progress. It takes as input a "guess" for
% the value of the permeabilities and gives as output the vector with the
% optimized permeabilities par_vec and the value of the objective function
% at the last iteration fval.

options = optimset('Display','final','PlotFcns',@optimplotfval,'TolFun',1e-5,...
    'TolX',1e-5);
[par_vec,fval]=fminsearch(@objective_nd,par_guess,options);
disp(fval)
end