% code for plotting the sensitivity indices in bar plots
clear all
close all
%%
load('X.mat');
load('Yfull.mat');

% Parameter Labels
efast_var={'$P_{NKCC}$','$P_{KIR,b}$','$P_{KIR,a}$',...
    '$P_{CFTR}$','$P_{N2B}$','$P_{N3B}$','$P_{AE,a}$','$P_{AE,b}$','$P_0^{TJ}$',...
    '$P_1^{TJ}$','$P_2^{TJ}$','$P_3^{TJ}$','dummy'};
% 13 permeabilities
cols=[1,2,3,4,5,6,7,8,9,10,11,12,13];

% numbering of the outputs
% 1'cc0',2'cc1',3'cc2',4'cc3',5'vc',6'va',7'xx',8'J0bc',9'J1bc',10'J2bc',
%    11'J3bc',12'J0ca',13'J1ca',14'J2ca',15'J3ca',16'Jtj',17'J1tj',18'J2tj',19'J3tj';

%%
% Sti plots 
% k refers to the outputs (see numbering above)
% choose the outputs you want and plot their sensitivity indices in a
% histogram
%%
k=1:6;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:13,[Sti(cols,k)],'Linewidth',2)  % Sti
set(gca,'xticklabel',efast_var(cols),'Linewidth',2)
legend('$c_0^c$','$c_1^c$','$c_2^c$','$c_3^c$','$V^c$','$V^a$')
% title('Cell concentrations and potential')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'STI_conc_pot_cell')
% exportgraphics(gcf,'STI_conc_pot_cell.png','Resolution',200)

%%
k=8:11;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:13,[Sti(cols,k)],'Linewidth',2)
set(gca,'xticklabel',efast_var(cols),'Linewidth',2)
legend('$\mathcal{J}_0^{bc}$','$\mathcal{J}_1^{bc}$','$\mathcal{J}_2^{bc}$','$\mathcal{J}_3^{bc}$')
%title('Basal to cell fluxes')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'STI_bc_fluxes_cell')
% exportgraphics(gcf,'STI_bc_fluxes_cell.png','Resolution',200)

%%
k=12:15;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:13,[Sti(cols,k)],'Linewidth',2)
set(gca,'xticklabel',efast_var(cols),'Linewidth',2)
legend('$\mathcal{J}_0^{ca}$','$\mathcal{J}_1^{ca}$','$\mathcal{J}_2^{ca}$','$\mathcal{J}_3^{ca}$')
%title('Cell to apical fluxes')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'STI_ca_fluxes_cell')
% exportgraphics(gcf,'STI_ca_fluxes_cell.png','Resolution',200)

%% 
k=16:19;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:13,[Sti(cols,k)],'Linewidth',2)
set(gca,'xticklabel',efast_var(cols),'Linewidth',2)
legend('$\mathcal{J}_0^{TJ}$','$\mathcal{J}_1^{TJ}$','$\mathcal{J}_2^{TJ}$','$\mathcal{J}_3^{TJ}$')
%title('Tight junction fluxes')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'STI_tj_fluxes_cell')
% exportgraphics(gcf,'STI_tj_fluxes_cell.png','Resolution',200)

%%
% Si plots 
% k refers to the outputs (see numbering above)

%%
k=1:6;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:13,[Si(cols,k)],'Linewidth',2)  % Sti
set(gca,'xticklabel',efast_var(cols),'Linewidth',2)
legend('$c_0^c$','$c_1^c$','$c_2^c$','$c_3^c$','$V_c$','$V_a$')
%title('Cell concentrations and potential')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'SI_conc_pot')

%% 
k=8:11;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:13,[Si(cols,k)],'Linewidth',2)
set(gca,'xticklabel',efast_var(cols),'Linewidth',2)
legend('$\mathcal{J}_0^{bc}$','$\mathcal{J}_1^{bc}$','$\mathcal{J}_2^{bc}$','$\mathcal{J}_3^{bc}$')
%title('Basal to cell fluxes')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'SI_bc')

%%
k=12:15;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:13,[Si(cols,k)],'Linewidth',2)
set(gca,'xticklabel',efast_var(cols),'Linewidth',2)
legend('$\mathcal{J}_0^{ca}$','$\mathcal{J}_1^{ca}$','$\mathcal{J}_2^{ca}$','$\mathcal{J}_3^{ca}$')
%title('Cell to apical fluxes')
set(gca,'Linewidth',2,'Fontsize',14)
% saveas(gcf,'SI_ca')

%%
k=16:19;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:13,[Si(cols,k)],'Linewidth',2)
set(gca,'xticklabel',efast_var(cols),'Linewidth',2)
legend('$\mathcal{J}_0^{tj}$','$\mathcal{J}_1^{tj}$','$\mathcal{J}_2^{tj}$','$\mathcal{J}_3^{tj}$')
%title('Tight junction fluxes')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'SI_tj')

%% check
% check the sum of the STi indices is bigger than 1
k=1:19;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);
sum_check_Sti=zeros(1,19);
for i=1:19
    sum_check_Sti(i)=sum(Sti(cols,i));
end
figure
plot(1:19,sum_check_Sti,'ko')
hold on
grid on
plot(1:19,ones(1,19),'r')
title('Sum total sensitivity indeces for each output') % should be >1
set(gca,'Fontsize',14)

% check the sum of the Si indices is smaller than 1
k=1:19;
sum_check_Si=zeros(1,19);
for i=1:19
    sum_check_Si(i)=sum(Si(cols,i));
end
figure
plot(1:19,sum_check_Si,'ko')
hold on
grid on
plot(1:19,ones(1,19),'r')
title('Sum first order sensitivity indeces for each output') % should be <1
set(gca,'Fontsize',14)
