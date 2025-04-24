
clear all 
close all

global  F R T Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx Ptj Atj P1tj P2tj P3tj...
 ca0 ca1 ca2 ca3  ...
 cb0 cb1 cb2 cb3 vb   ...
 L h Ux Q C...
 xmesh w...
 D0 D1 D2 D3 ...
 w_cheb_01 ...
 ref epsilon mu lambdaD chi bar_sigma0 H fac half_fac npts D Lp Ltj nit

% select values of Lp_vec to test
Lp=2e-12;
Ltj_mult=10;

sigma0_vec=linspace(-0.0243,0,20); % wall charge
Cm=1e-2;

%% EO present. Prepare matrices for the values to save:
%% SKIP

EO_wf_matrix=zeros(length(sigma0_vec),1); % transendothelial water flux
EO_osm_press_la_dim=zeros(length(sigma0_vec),1); % osmotic pressure
EO_mech_press_tj_dim=zeros(length(sigma0_vec),1); % mechanical pressure
EO_Va=zeros(length(sigma0_vec),1); % transendothelial potential
EO_Vc=zeros(length(sigma0_vec),1); % cellular potential

% cellular concentrations
EO_Cc0=zeros(length(sigma0_vec),1); 
EO_Cc1=zeros(length(sigma0_vec),1);
EO_Cc2=zeros(length(sigma0_vec),1);
EO_Cc3=zeros(length(sigma0_vec),1);

% transendothelial ion fluxes
EO_Na_flux=zeros(length(sigma0_vec),1); % J0tj+J0ca
EO_K_flux=zeros(length(sigma0_vec),1);
EO_Cl_flux=zeros(length(sigma0_vec),1);
EO_HCO3_flux=zeros(length(sigma0_vec),1);

% cleft concentrations at the tight junction
EO_Cl0=zeros(length(sigma0_vec),1); 
EO_Cl1=zeros(length(sigma0_vec),1);
EO_Cl2=zeros(length(sigma0_vec),1);
EO_Cl3=zeros(length(sigma0_vec),1);


%% compute the solution for all the selected Lp values
%% SKIP 

for j=1:length(sigma0_vec)
        sigma0=sigma0_vec(j); 
        Ltj=Lp*Ltj_mult;
        % run the model
        main_coupling
        % save the results
        EO_wf_matrix(j)=overall_water_flux;
        EO_osm_press_la_dim(j)= R*T*(ca0+ca1+ca2+ca3-cl0(end)-cl1(end)-cl2(end)-cl3(end))*C;
        EO_mech_press_tj_dim(j)=y9(end)*mu*Ux*L/h^2;
        EO_Cc0(j)=Cc0;
        EO_Cc1(j)=Cc1;
        EO_Cc2(j)=Cc2;
        EO_Cc3(j)=Cc3;
        EO_Vc(j)=Vc;
        EO_Va(j)=Va;
        EO_Na_flux(j)=J0ca+Jtj; 
        EO_K_flux(j)=J1ca+J1tj;
        EO_Cl_flux(j)=J2ca+J2tj;
        EO_HCO3_flux(j)=J3ca+J3tj;
        EO_Cl0(j)=cl0(end)*C; 
        EO_Cl1(j)=cl1(end)*C;
        EO_Cl2(j)=cl2(end)*C;
        EO_Cl3(j)=cl3(end)*C;

end
save('EO_wf_matrix.mat','EO_wf_matrix');
save('EO_mech_press_tj_dim','EO_mech_press_tj_dim');
save('EO_osm_press_la_dim','EO_osm_press_la_dim');
save('EO_Cc0.mat','EO_Cc0')
save('EO_Cc1.mat','EO_Cc1')
save('EO_Cc2.mat','EO_Cc2')
save('EO_Cc3.mat','EO_Cc3')
save('EO_Vc.mat','EO_Vc')
save('EO_Va.mat','EO_Va')
save('EO_Na_flux.mat','EO_Na_flux')
save('EO_K_flux.mat','EO_K_flux')
save('EO_Cl_flux.mat','EO_Cl_flux')
save('EO_HCO3_flux.mat','EO_HCO3_flux')
save('EO_Cl0.mat','EO_Cl0')
save('EO_Cl1.mat','EO_Cl1')
save('EO_Cl2.mat','EO_Cl2')
save('EO_Cl3.mat','EO_Cl3')


%% run from here for plots if you already have the matrices

%% EO absent, only local osmosis

clear all

global  F R T Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx Ptj Atj P1tj P2tj P3tj...
 ca0 ca1 ca2 ca3  ...
 cb0 cb1 cb2 cb3 vb   ...
 L h Ux Q C...
 xmesh w...
 D0 D1 D2 D3 ...
 w_cheb_01 ...
 ref epsilon mu lambdaD chi bar_sigma0 H fac half_fac npts D Lp Ltj

Lp=2e-12;
Ltj_mult=10;

sigma0=0;
Cm=0;


% compute the model solution for all the values of Lp
        Ltj=Lp*Ltj_mult;
        
        % solve the model
        main_coupling
        % save the values
    LO_flux=overall_water_flux;
       
%% load matrices for the plots 
% clear all 

% EO

load('EO_mech_press_tj_dim.mat');
load('EO_osm_press_la_dim.mat');
load('EO_wf_matrix.mat');
load('EO_Cc0.mat');
load('EO_Cc1.mat');
load('EO_Cc2.mat');
load('EO_Cc3.mat');
load('EO_Vc.mat');
load('EO_Va.mat');
load('EO_Na_flux.mat');
load('EO_Cl_flux.mat');
load('EO_K_flux.mat');
load('EO_HCO3_flux.mat');
load('EO_Cl0.mat');
load('EO_Cl1.mat');
load('EO_Cl2.mat');
load('EO_Cl3.mat');

% %% water flux across the endothelium with and without EO
% sigma0_vec=linspace(-0.0243,0,20);
% figure('WindowState','maximized')
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaultLegendInterpreter','latex');
% % loglog(Lp_vec,noEO_wf_matrix,'k','LineWidth',3)
% % hold on
% % loglog(Lp_vec,EO_wf_matrix-noEO_wf_matrix,'r','LineWidth',3)
% % semilogy(sigma0_vec(1:end-1),EO_wf_matrix(1:end-1)-EO_wf_matrix(end),'r','LineWidth',3)
% % plot(sigma0_vec,EO_wf_matrix,'r','LineWidth',3)
% plot(sigma0_vec,EO_wf_matrix-EO_wf_matrix(end),'r','LineWidth',3)
% 
% xlabel('$\sigma_0$ (C/m$^2$)','Interpreter','latex')
% xlim([sigma0_vec(1) sigma0_vec(end)])
% ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
% xline(-0.0048,'--k','LineWidth',3)
% %leg=legend('Local Osmosis','Electro-osmosis');
% %leg.Location="southeast";
% % title('Water flux across the endothelium with and without EO')
% set(gca,'LineWidth',3, 'FontSize',28)
% % saveas(gcf,'Water flux')
% % exportgraphics(gcf,'EOvsLO.png','Resolution',200)
% 
% 
% perc_diff=(EO_wf_matrix(1)-EO_wf_matrix(end))/EO_wf_matrix(end)*100;

%%
sigma0_vec=linspace(-0.0243,0,20);
figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
% loglog(Lp_vec,noEO_wf_matrix,'k','LineWidth',3)
% hold on
% loglog(Lp_vec,EO_wf_matrix-noEO_wf_matrix,'r','LineWidth',3)
% semilogy(sigma0_vec(1:end-1),EO_wf_matrix(1:end-1)-EO_wf_matrix(end),'r','LineWidth',3)
% plot(sigma0_vec,EO_wf_matrix,'r','LineWidth',3)
plot(sigma0_vec,EO_wf_matrix-LO_flux,'r','LineWidth',3)


xlabel('$\sigma_0$ (C/m$^2$)','Interpreter','latex')
xlim([sigma0_vec(1) sigma0_vec(end)])
ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
xline(-0.0048,'--k','LineWidth',3)
%leg=legend('Local Osmosis','Electro-osmosis');
%leg.Location="southeast";
% title('Water flux across the endothelium with and without EO')
set(gca,'LineWidth',3, 'FontSize',28)
% saveas(gcf,'Water flux')
exportgraphics(gcf,'EO_flux_vs_sigma0.png','Resolution',200)


perc_diff=(EO_wf_matrix(1)-LO_flux)/LO_flux*100;
