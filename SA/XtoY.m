% this script computes the model output for each of the samples of the
% parameter space in the matrix X.

clear;
close all;

global C ca0 ca1 ca2 ca3 cb0 cb1 cb2 cb3 vb npts

load('X.mat') % call X file

Parameter_settings_EFAST; % call parameter file

Y(NS,length(y0),length(pmin),NR)=0;  % pre-allocation

for i=1:k
for Ll=1:NR
        for run_num=1:NS
            [i run_num Ll] % keeps track of [parameter run NR]
            % ODE system file
            [cc0,cc1,cc2,cc3,vc,va,xx,J0tot,J1tot,J2tot,J3tot,Jtj,J1tj,J2tj,J3tj,overall_water_flux,u_end_hex]=model_main(X(:,:,i,Ll),run_num);
            % It saves only the desired output
            Y(run_num,:,i,Ll)=[cc0,cc1,cc2,cc3,vc,va,xx,J0tot,J1tot,J2tot,J3tot,Jtj,J1tj,J2tj,J3tj,overall_water_flux]'; % this is the output 
            save res.mat;
        end % run_num=1:NS
end % L=1:NR
end

save res.mat;
% exit
save('Yfull.mat','Y')
