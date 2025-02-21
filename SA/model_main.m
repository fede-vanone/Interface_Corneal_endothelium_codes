function [cc0,cc1,cc2,cc3,vc,va,xx,J0tot,J1tot,J2tot,J3tot,Jtj,J1tj,J2tj,J3tj,overall_water_flux,u_end_hex]=model_main(X,run_num)
% this function takes as input the matrix X (containing the samples of the 
% parameter space) and run_num (the number of the iteration in the loop of 
% XtoY.m) and runs main_coupling.m with the corresponding permeability values.
global  F R T Lp Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx Ptj Atj P1tj P2tj P3tj...
 ca0 ca1 ca2 ca3 Aa ...
 cb0 cb1 cb2 cb3 vb Ab  ...
 L h Ux Q C...
 xmesh w...
 D0 D1 D2 D3 ...
 Ltj w_cheb_01 ...
 ref epsilon mu lambdaD sigma0 Cm chi bar_sigma0 H fac half_fac npts D
Parameter_settings_EFAST;

Pnkcc = X(run_num,1);
Pkirb = X(run_num,2);
Pkira = X(run_num,3);
Pcftr = X(run_num,4);
Pn2b = X(run_num,5);
Pn3b = X(run_num,6);
Ptj=X(run_num,7);
P1tj = X(run_num,8);
P2tj= X(run_num,9);
P3tj= X(run_num,10);
dummy = X(run_num,11);

try
 main_coupling; % script to run with parameters above 

catch
    cc0=NaN;cc1=NaN;cc2=NaN;cc3=NaN;vc=NaN;va=NaN;xx=NaN;J0tot=NaN;J1tot=NaN;...
        J2tot=NaN;J3tot=NaN;Jtj=NaN;...
        J1tj=NaN;J2tj=NaN;J3tj=NaN;overall_water_flux=NaN;u_end_hex=NaN;

end