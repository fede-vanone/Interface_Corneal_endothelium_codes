% code for plotting the sensitivity indices in bar plots
clear all
close all
%%
load('X.mat');
load('Yfull.mat');

% Parameter Labels
efast_var={'$P_{NKCC}$','$P_{KIR,b}$','$P_{KIR,a}$',...
    '$P_{CFTR}$','$P_{N2B}$','$P_{N3B}$','$P_0^{TJ}$',...
    '$P_1^{TJ}$','$P_2^{TJ}$','$P_3^{TJ}$','dummy'};
% 11 permeabilities
cols=[1,2,3,4,5,6,7,8,9,10,11];

% numbering of the outputs
% 1 cc0,2 cc1,3 cc2,4 cc3,5 vc,6 va,7 xx,8 J0tot,9 J1tot,10 J2tot,11 J3tot,
% 12 Jtj,13 J1tj,14 J2tj,15 J3tj,16 overall_water_flux,17 u_end_hex

%%
% Sti plots 
% k refers to the outputs (see numbering above)
% choose the outputs you want and plot their sensitivity indices in a
% histogram

%%
k=1:6;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

% figure('units','normalized','outerposition',[0 0 1 1])
figure('WindowState','maximized');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Sti(cols,k)],'LineWidth',2)  % Sti
set(gca,'xticklabel',efast_var(cols),'LineWidth',2)
legend('$c_0^c$','$c_1^c$','$c_2^c$','$c_3^c$','$V_c$','$V_a$')
% title('Cell concentrations and potential')
set(gca,'Linewidth',2,'Fontsize',24,'FontWeight','bold')
% saveas(gcf,'STI_conc_pot')
% exportgraphics(gcf,'STI_conc_pot_res200.png','Resolution',200)

%%
k=8:11;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Sti(cols,k)],'LineWidth',2)
set(gca,'xticklabel',efast_var(cols),'LineWidth',2,'FontWeight','bold')
legend('$\mathcal{J}_0$','$\mathcal{J}_1$','$\mathcal{J}_2$','$\mathcal{J}_3$')
% title('Overall ion fluxes across the endothelium')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'STI_tot_fluxes')
% exportgraphics(gcf,'STI_tot_fluxes_res200.png','Resolution',200)

%%
k=12:15;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Sti(cols,k)],'LineWidth',2)
set(gca,'xticklabel',efast_var(cols),'LineWidth',2)
legend('$\mathcal{J}_0^{TJ}$','$\mathcal{J}_1^{TJ}$','$\mathcal{J}_2^{TJ}$','$\mathcal{J}_3^{TJ}$')
% title('Tight junction ion fluxes')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'STI_tj_fluxes')
% exportgraphics(gcf,'STI_tj_fluxes_res200.png','Resolution',200)

%%
k=16;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Sti(cols,k)],'LineWidth',2)
set(gca,'xticklabel',efast_var(cols),'LineWidth',2)
legend('$\bar{Q}$')
% title('Overall water flux across the endothelium')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'STI_water_flux')
% exportgraphics(gcf,'STI_water_flux_res200.png','Resolution',200)

%%
k=17;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Sti(cols,k)],'LineWidth',2)
set(gca,'xticklabel',efast_var(cols),'LineWidth',2)
legend('$\bar{u}$')
% title('Overall water flux across the endothelium')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'STI_water_flux_hexagon')
% exportgraphics(gcf,'STI_water_flux_hex_res200.png','Resolution',200)

%%
% Si plots 
% k refers to the outputs (see numbering above)

%%
k=1:4;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Sti(cols,k)],'LineWidth',2) % 
set(gca,'xticklabel',efast_var(cols),'LineWidth',2)
legend('$c_0^c$','$c_1^c$','$c_2^c$','$c_3^c$')
% title('Total sensitivity indices - Cell concentrations')
set(gca,'Linewidth',2,'Fontsize',24)

%%
k=5:6;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure('WindowState','maximized');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Sti(cols,k)],'LineWidth',2) % 
set(gca,'xticklabel',efast_var(cols),'LineWidth',2)
legend('$V^c$','$V^a$')
% title('Total sensitivity indices - Potential')
set(gca,'Linewidth',2,'Fontsize',24)
% saveas(gcf,'STI_pot')
% exportgraphics(gcf,'STI_pot_res200.png','Resolution',200)

%%
k=8:11;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Si(cols,k)],'LineWidth',2)
set(gca,'xticklabel',efast_var(cols),'LineWidth',2,'FontWeight','bold')
legend('$\mathbf{\mathcal{J}_0}$','$\mathbf{\mathcal{J}_1}$','$\mathbf{\mathcal{J}_2}$',...
    '$\mathbf{\mathcal{J}_3}$')
title('Overall ion fluxes across the endothelium')
set(gca,'Linewidth',2,'Fontsize',24,'FontWeight','bold')

%%
k=12:15;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Si(cols,k)],'LineWidth',2)
set(gca,'xticklabel',efast_var(cols),'LineWidth',2,'FontWeight','bold')
legend('$\mathbf{\mathcal{J}_0^{TJ}}$','$\mathbf{\mathcal{J}_1^{TJ}}$',...
    '$\mathbf{\mathcal{J}_2^{TJ}}$','$\mathbf{\mathcal{J}_3^{TJ}}$')
title('Tight junction ion fluxes')
set(gca,'Linewidth',2,'Fontsize',24,'FontWeight','bold')

%%
k=16;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
bar(1:11,[Si(cols,k)],'LineWidth',2)
set(gca,'xticklabel',efast_var(cols))
legend('$\mathbf{Q^{te}}$')
title('Overall water flux across the endothelium')
set(gca,'Linewidth',2,'Fontsize',24,'FontWeight','bold')

%% check
% check the sum of the STi indices is bigger than 1

k=1:17;
[Si,Sti,rangeSi,rangeSti] = efast_sd_2(Y,OMi,MI,k);
sum_check_Sti=zeros(1,17);
for i=1:17
    sum_check_Sti(i)=sum(Sti(cols,i));
end
figure
plot(1:17,sum_check_Sti,'ko')
hold on
grid on
plot(1:17,ones(1,17),'r')
title('Sum total sensitivity indeces for each output') % should be >1
set(gca,'Fontsize',14)
saveas(gcf,'checkSTI')

% check the sum of the Si indices is smaller than 1
k=1:17;
sum_check_Si=zeros(1,17);
for i=1:17
    sum_check_Si(i)=sum(Si(cols,i));
end
figure
plot(1:17,sum_check_Si,'ko')
hold on
grid on
plot(1:17,ones(1,17),'r')
title('Sum first order sensitivity indeces for each output') % should be <1
set(gca,'Fontsize',14)
saveas(gcf,'checkSI')


