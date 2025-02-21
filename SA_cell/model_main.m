function [cc0,cc1,cc2,cc3,vc,va,xx,J0bc,J1bc,J2bc,J3bc,J0ca,J1ca,J2ca,J3ca,Jtj,J1tj,J2tj,J3tj]=model_main(X,run_num)
% this function takes as input the matrix X (containing the samples of the 
% parameter space) and run_num (the number of the iteration in the loop of 
% XtoY.m) and runs cell_full.m with the corresponding permeability values. 
global C ca0 ca1 ca2 ca3 cb0 cb1 cb2 cb3 vb npts
Parameter_settings_EFAST;

Pnkcc = X(run_num,1);
Pkirb = X(run_num,2);
Pkira = X(run_num,3);
Pcftr = X(run_num,4);
Pn2b = X(run_num,5);
Pn3b = X(run_num,6);
Paea=X(run_num,7);
Pae = X(run_num,8);
Ptj=X(run_num,9);
P1tj=X(run_num,10);
P2tj=X(run_num,11);
P3tj=X(run_num,12);
dummy = X(run_num,13);


try
cell_full; % script to run with parameters above 
catch

cc0=NaN; cc1=NaN; cc2=NaN; cc3=NaN; vc=NaN; va=NaN; xx=NaN; J0bc=NaN; J1bc=NaN; ...
    J2bc=NaN; J3bc=NaN; J0ca=NaN; J1ca=NaN; J2ca=NaN; J3ca=NaN; ...
    Jtj=NaN; J1tj=NaN; J2tj=NaN; J3tj=NaN;
end