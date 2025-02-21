
% EO
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

Lp=2e-12;
Ltj_mult_vec=10.^linspace(0,3,20);

sigma0=-0.0048; % wall charge
Cm=1e-2;

%%
EO_wf_hexagon=zeros(length(Ltj_mult_vec),1);
EO_wf_matrix=zeros(length(Ltj_mult_vec),1);

EO_osm_press_la_dim=zeros(length(Ltj_mult_vec),1);
EO_mech_press_tj_dim=zeros(length(Ltj_mult_vec),1);

EO_Cc0=zeros(length(Ltj_mult_vec),1);
EO_Cc1=zeros(length(Ltj_mult_vec),1);
EO_Cc2=zeros(length(Ltj_mult_vec),1);
EO_Cc3=zeros(length(Ltj_mult_vec),1);

EO_Va=zeros(length(Ltj_mult_vec),1);
EO_Vc=zeros(length(Ltj_mult_vec),1);

EO_Na_flux=zeros(length(Ltj_mult_vec),1); % J0tj+J0ca
EO_K_flux=zeros(length(Ltj_mult_vec),1);
EO_Cl_flux=zeros(length(Ltj_mult_vec),1);
EO_HCO3_flux=zeros(length(Ltj_mult_vec),1);

EO_Cl0=zeros(length(Ltj_mult_vec),1); % cleft concentration at the tight junction
EO_Cl1=zeros(length(Ltj_mult_vec),1);
EO_Cl2=zeros(length(Ltj_mult_vec),1);
EO_Cl3=zeros(length(Ltj_mult_vec),1);




%%
for j=1:length(Ltj_mult_vec)
    % clear Lp
        Ltj_mult=Ltj_mult_vec(j);
        Ltj=Lp*Ltj_mult;
        
        main_coupling
        EO_wf_hexagon(j)=u_end_hex;
        EO_wf_matrix(j)=overall_water_flux;
        EO_osm_press_la_dim(j)= R*T*(ca0+ca1+ca2+ca3-cl0(end)-cl1(end)-cl2(end)-cl3(end))*C;
        EO_mech_press_tj_dim(j)=y9(end)*mu*Ux*L/h^2;
        EO_Cc0(j)=Cc0;
        EO_Cc1(j)=Cc1;
        EO_Cc2(j)=Cc2;
        EO_Cc3(j)=Cc3;
        EO_Vc(j)=Vc;
        EO_Va(j)=Va;
        EO_Na_flux(j)=J0ca+Jtj; % J0tj+J0ca
        EO_K_flux(j)=J1ca+J1tj;
        EO_Cl_flux(j)=J2ca+J2tj;
        EO_HCO3_flux(j)=J3ca+J3tj;
        EO_Cl0(j)=cl0(end)*C; % cleft concentration at the tight junction
        EO_Cl1(j)=cl1(end)*C;
        EO_Cl2(j)=cl2(end)*C;
        EO_Cl3(j)=cl3(end)*C;

end
save('EO_wf_hex.mat','EO_wf_hexagon')
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

%%

% noEO
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
Ltj_mult_vec=10.^linspace(0,3,20);

sigma0=0; % wall charge
Cm=0;

noEO_wf_hexagon=zeros(length(Ltj_mult_vec),1);
noEO_wf_matrix=zeros(length(Ltj_mult_vec),1);
noEO_osm_press_la_dim=zeros(length(Ltj_mult_vec),1);
noEO_mech_press_tj_dim=zeros(length(Ltj_mult_vec),1);

noEO_Cc0=zeros(length(Ltj_mult_vec),1);
noEO_Cc1=zeros(length(Ltj_mult_vec),1);
noEO_Cc2=zeros(length(Ltj_mult_vec),1);
noEO_Cc3=zeros(length(Ltj_mult_vec),1);

noEO_Va=zeros(length(Ltj_mult_vec),1);
noEO_Vc=zeros(length(Ltj_mult_vec),1);

noEO_Na_flux=zeros(length(Ltj_mult_vec),1); % J0tj+J0ca
noEO_K_flux=zeros(length(Ltj_mult_vec),1);
noEO_Cl_flux=zeros(length(Ltj_mult_vec),1);
noEO_HCO3_flux=zeros(length(Ltj_mult_vec),1);

noEO_Cl0=zeros(length(Ltj_mult_vec),1); % cleft concentration at the tight junction
noEO_Cl1=zeros(length(Ltj_mult_vec),1);
noEO_Cl2=zeros(length(Ltj_mult_vec),1);
noEO_Cl3=zeros(length(Ltj_mult_vec),1);

% Lp_vec=[2e-13 2e-12 1e-11 2e-11];
% Ltj_mult=[0 1 10 100 1000];

for j=1:length(Ltj_mult_vec)
    % clear Lp
        Ltj_mult=Ltj_mult_vec(j);
        Ltj=Lp*Ltj_mult;
        
        main_coupling
        noEO_wf_hexagon(j)=u_end_hex;
        noEO_wf_matrix(j)=overall_water_flux;
        noEO_osm_press_la_dim(j)= R*T*(ca0+ca1+ca2+ca3-cl0(end)-cl1(end)-cl2(end)-cl3(end))*C;
        noEO_mech_press_tj_dim(j)=y9(end)*mu*Ux*L/h^2;
        noEO_Cc0(j)=Cc0;
        noEO_Cc1(j)=Cc1;
        noEO_Cc2(j)=Cc2;
        noEO_Cc3(j)=Cc3;
        noEO_Vc(j)=Vc;
        noEO_Va(j)=Va;
        noEO_Na_flux(j)=J0ca+Jtj; % J0tj+J0ca
        noEO_K_flux(j)=J1ca+J1tj;
        noEO_Cl_flux(j)=J2ca+J2tj;
        noEO_HCO3_flux(j)=J3ca+J3tj;
        noEO_Cl0(j)=cl0(end)*C; % cleft concentration at the tight junction
        noEO_Cl1(j)=cl1(end)*C;
        noEO_Cl2(j)=cl2(end)*C;
        noEO_Cl3(j)=cl3(end)*C;

end

save('noEO_wf_end.mat','noEO_wf_hexagon')
save('noEO_wf_matrix.mat','noEO_wf_matrix');
save('noEO_mech_press_tj_dim.mat','noEO_mech_press_tj_dim');
save('noEO_osm_press_la_dim.mat','noEO_osm_press_la_dim');
save('noEO_Cc0.mat','noEO_Cc0')
save('noEO_Cc1.mat','noEO_Cc1')
save('noEO_Cc2.mat','noEO_Cc2')
save('noEO_Cc3.mat','noEO_Cc3')
save('noEO_Vc.mat','noEO_Vc')
save('noEO_Va.mat','noEO_Va')
save('noEO_Na_flux.mat','noEO_Na_flux')
save('noEO_K_flux.mat','noEO_K_flux')
save('noEO_Cl_flux.mat','noEO_Cl_flux')
save('noEO_HCO3_flux.mat','noEO_HCO3_flux')
save('noEO_Cl0.mat','noEO_Cl0')
save('noEO_Cl1.mat','noEO_Cl1')
save('noEO_Cl2.mat','noEO_Cl2')
save('noEO_Cl3.mat','noEO_Cl3')

%% run from here for plots if you already have the matrices
%% load matrices for the plots 
clear all 

% EO
Lp=2e-12;
Ltj_mult_vec=10.^linspace(0,3,20);

load('EO_wf_hex.mat')
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

% noEO

load("noEO_wf_end.mat") % hexagon
load('noEO_mech_press_tj_dim.mat');
load('noEO_osm_press_la_dim.mat');
load('noEO_wf_matrix.mat');
load('noEO_Cc0.mat');
load('noEO_Cc1.mat');
load('noEO_Cc2.mat');
load('noEO_Cc3.mat');
load('noEO_Vc.mat');
load('noEO_Va.mat');
load('noEO_Na_flux.mat');
load('noEO_Cl_flux.mat');
load('noEO_K_flux.mat');
load('noEO_HCO3_flux.mat');
load('noEO_Cl0.mat');
load('noEO_Cl1.mat');
load('noEO_Cl2.mat');
load('noEO_Cl3.mat');

%% water flux across the endothelium with and without EO
% close(1)

figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
semilogy(Ltj_mult_vec,noEO_wf_matrix,'k','LineWidth',3)
hold on
semilogy(Ltj_mult_vec,EO_wf_matrix-noEO_wf_matrix,'r','LineWidth',3)
% xlabel('$L_p$ (m/Pa/s)','Interpreter','latex')
xlim([Ltj_mult_vec(1) Ltj_mult_vec(end)])
ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
% xline(2e-12,'--k','LineWidth',3)
leg=legend('Local Osmosis','Electro-osmosis');
leg.Location="southeast";
% title('Water flux across the endothelium with and without EO')

set(gca,'LineWidth',3, 'FontSize',24)
% saveas(gcf,'Water flux Ltj ')
% exportgraphics(gcf,'EOvsLO_Ltj.png','Resolution',200)


%% log plot Water flux EO vs noEO
figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
loglog(Ltj_mult_vec,noEO_wf_matrix,'k','LineWidth',3)
hold on
loglog(Ltj_mult_vec,EO_wf_matrix-noEO_wf_matrix,'r','LineWidth',3)
xlabel('$L_{TJ}/L_p$','Interpreter','latex')
xlim([Ltj_mult_vec(1) Ltj_mult_vec(end)])
ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
xline(10,'--k','LineWidth',3)
leg=legend('Local Osmosis','Electro-osmosis');
leg.Location="southeast";
% title('Water flux across the endothelium with and without EO')

set(gca,'LineWidth',3, 'FontSize',28)
saveas(gcf,'Water flux Ltj log')
exportgraphics(gcf,'EOvsLO_Ltj_log.png','Resolution',200)

%% mechanical pressure at the tj with and without EO
% close(2)
figure(2)
plot(Ltj_mult_vec,EO_mech_press_tj_dim,'g--','LineWidth',2)
hold on
plot(Ltj_mult_vec,noEO_mech_press_tj_dim,'g','LineWidth',2)
legend('EO','LO')
ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
ylabel('Pascal')
% title('Mechanical pressure at the tj with and without EO')
set(gca,'linewidth',2, 'FontSize',16, 'FontWeight','bold')
% saveas(gcf,'Mechanical pressure')

%% mechanical and osmotic pressure at the tight junction with EO
% close(3)
figure(3)
plot(Ltj_mult_vec,EO_mech_press_tj_dim,'g--','LineWidth',2)
hold on
plot(Ltj_mult_vec,EO_osm_press_la_dim,'g','LineWidth',2)
legend('Mechanical','Osmotic')
ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
ylabel('Pressure (Pa)')

% title('Mechanical and osmotic pressure at the tight junction with EO')
set(gca,'linewidth',2, 'FontSize',16, 'FontWeight','bold')
% saveas(gcf,'Osm vs mech EO')

%% difference in water flux OE-noEO
%close(4)
figure(4)
plot(Ltj_mult_vec,EO_wf_matrix-noEO_wf_matrix,'g--','LineWidth',2)
hold on
ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
ylabel('m/s')
% title('Water flux difference across the endothelium EO-noEO')
set(gca,'linewidth',2, 'FontSize',16, 'FontWeight','bold')
% saveas(gcf,'Water flux difference Eo vs LO')

%% percentage difference water flux with and without EO
% close(5)
figure(5)
plot(Ltj_mult_vec,(EO_wf_matrix-noEO_wf_matrix)./EO_wf_matrix*100,'g--','LineWidth',2)
ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
% ylabel('m/s')
% title('Water flux difference across the endothelium EO-noEO')
set(gca,'linewidth',2, 'FontSize',16, 'FontWeight','bold')
% saveas(gcf,'Perc water flux difference EO LO')

%% mechanical and osmotic pressure at the tight junction without EO
% close(6)
figure(6)
plot(Ltj_mult_vec,noEO_mech_press_tj_dim,'g--','LineWidth',2)
hold on
plot(Ltj_mult_vec,noEO_osm_press_la_dim,'g','LineWidth',2)
legend('Mechanical','Osmotic')
ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
ylabel('Pressure (Pa)')
% title('Mechanical and osmotic pressure at the tight junction without EO')
set(gca,'linewidth',2, 'FontSize',16, 'FontWeight','bold')
% saveas(gcf,'OSM vs mech LO')

%% osmotic pressure at the tj with and without EO
% close(7)
figure(7)
plot(Ltj_mult_vec,EO_osm_press_la_dim,'g--','LineWidth',2)
hold on
plot(Ltj_mult_vec,noEO_osm_press_la_dim,'g','LineWidth',2)
legend('EO','LO')
ylabel('$\bar{Q}$ (m/s)','Interpreter','latex')
ylabel('Pressure (Pa)')
% title('Osmotic pressure at the tj with and without EO')
set(gca,'linewidth',2, 'FontSize',16, 'FontWeight','bold')
% saveas(gcf,'Osmotic pressure')

%% difference in water flux considering or not the hexagonal shape of cells
figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
plot(EO_wf_matrix,EO_wf_hexagon,'LineWidth',3)
xlabel('$\bar{Q}$ (m/s)','Interpreter','latex')
xlim([EO_wf_matrix(1) EO_wf_matrix(end)])
ylabel('$\bar{u}$ (m/s)','Interpreter','latex')

set(gca,'LineWidth',3, 'FontSize',24)
% saveas(gcf,'Water flux Ltj ')
% print(gcf,'Qvsu_hexagon.png','-dpng', '-r200')


sum(EO_osm_press_la_dim-EO_mech_press_tj_dim)/length(EO_osm_press_la_dim);
sum(EO_mech_press_tj_dim-noEO_mech_press_tj_dim)/length(EO_mech_press_tj_dim);


