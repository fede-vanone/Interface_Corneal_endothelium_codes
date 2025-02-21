
% this code runs the coupled model cell-cleft

clear all
close all

global  F R T Lp Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx Ptj Atj P1tj P2tj P3tj...
 ca0 ca1 ca2 ca3 Aa ...
 cb0 cb1 cb2 cb3 vb Ab  ...
 L h Ux Q C...
 xmesh w...
 D0 D1 D2 D3 ...
 Ltj w_cheb_01 ...
 ref epsilon mu lambdaD sigma0 Cm chi bar_sigma0 H fac half_fac npts D

%% Parameters

[F,R,T,h,L,Ux,Q,C,zx,Lp,Ltj,...
    Atj,w,D0,D1,D2,D3,cb0,cb1,cb2,cb3,vb,ca0,ca1,ca2,ca3,...
    toll,nit,npts,epsilon,mu,lambdaD,sigma0,Cm]=pars5();

ref=1e-6;

% Peclet number 
Pe=L*Ux/D0;
% fprintf('Peclet number=%g\n',Pe)

%% Permeabilities

% set1 
Ppump= 4.197716594624170e-06; 
Pnkcc= 3.519978970519485e-06; 
Pkirb= 2.074799294625940e-07;
Pkira= 1.544574301723440e-07;  
Pcftr= 1.743852885136938e-07;  
Pn2b= 4.661932732935699e-06; 
Pn3b= 2.972204545461175e-06;  
Paea= 5.318806049015642e-07;  
Pae= 4.547240291363707e-07; 
Ptj= 0.002438293282052; 
P1tj= 0.003077082580040; 
P2tj= 0.001131927484887; 
P3tj= 0.001002449976630; 

% % set2
% Ppump= 4.197716594624170e-06; 
% Pnkcc= 2.903345431655278e-06; 
% Pkirb= 4.158586873235697e-08; 
% Pkira= 5.259824354117620e-07; 
% Pcftr= 1.460988750783429e-07; 
% Pn2b= 1.989694826595409e-06;
% Pn3b= 6.007145420019043e-07; 
% Paea= 2.266921537667212e-06; 
% Pae= 9.514330239191147e-08;
% Ptj= 0.007203826257399; 
% P1tj= 6.186441835851560e-04;
% P2tj=0.001649264908772; 
% P3tj= 0.004084852835289; 

%% initialization of errors

% error on variables in the cell
err_cell_norm=zeros(1,nit);
% error on variables in the cleft
err_cleft_norm=zeros(1,nit);
% error in the variables in the cell + cleft
err_cell_cleft=zeros(1,nit);
% bvp5c error
err_bvp=zeros(1,nit);
% residuals of cell equations
res_cell=zeros(1,nit);

%% mesh for the endothelium ODEs solution

[w_cheb,x_cheb]=clenshaw_curtis_p(npts,1); % x_cheb in [1,-1]
xmesh=flip(x_cheb+1)'/2; % flip [-1,1], shift to [0,2] and rescale to [0,1]
mapping_01=1/2;
w_cheb_01=flip(w_cheb')*mapping_01; % weights for integration formula in 01
[x_chebdif, DM] = chebdif(npts, 1); 
% df/dx = D*x, x column vector
D=flip(flip(DM/mapping_01,1),2); % derivatives in [0,1] (flip rows and columns)

%% initial lateral concentrations and potential input

% initialization of the cleft variables (concentration and potential)
% for the first iteration, the cleft is homogeneous (as the basal region)

cl0=cb0*ones(1,npts);
dcl0=zeros(1,npts);
cl1=cb1*ones(1,npts);
dcl1=zeros(1,npts);
cl2=cb2*ones(1,npts);
dcl2=zeros(1,npts);
cl3=ones(1,npts)*cb3;
dcl3=zeros(1,npts);
vl=ones(1,npts)*vb;
dvl=zeros(1,npts);
p=ones(1,npts)*0;
dp=ones(1,npts)*0;

%% start with one cell solution with the lateral variables initialized as in the basal
[cc0_start,cc1_start,cc2_start,cc3_start,vc_start,va_start,xx_start,J0bc,J1bc,J2bc,J3bc,J0ca,J1ca,J2ca,J3ca]=...
    cell_full(npts,cl0,cl1,cl2,cl3,vl,p);

cc0=cc0_start; cc1=cc1_start; cc2=cc2_start; cc3=cc3_start; xx=xx_start; ...
    va=va_start; vc=vc_start;

%% iteration 

for n=1:nit
   fprintf('Iteration %g \n',n) 

% first solve for the cleft
[bvperr_n,y1_n,y2_n,y3_n,y4_n,y5_n,y6_n,y7_n,y8_n,y9_n,y10_n] = ...
    solve_end_cleft(npts,xx,cc0,cc1,cc2,cc3,va,vc);

% cl3 from electroneutrality
cl3_n=y1_n+y3_n-y5_n;

cl0_n=y1_n; % Cl0
dcl0_n=y2_n; % dCl0
cl1_n=y3_n; % Cl1
dcl1_n=y4_n; % dCl1
cl2_n=y5_n; % Cl2
dcl2_n=y6_n; % dCl2
vl_n=y7_n; % V
dvl_n=y8_n; % dV
p_n=y9_n; % p
dp_n=y10_n; % dp

% then solve for the cell
[cc0_n, cc1_n,cc2_n, cc3_n, vc_n, va_n, xx_n,J0bc_n,J1bc_n,J2bc_n,J3bc_n,J0ca_n,J1ca_n,J2ca_n,J3ca_n]...
    = cell_full(npts,cl0_n,cl1_n,cl2_n,cl3_n,vl_n,p_n);

% compute the errors
% relative cell error 
err_cell_vec=[abs(cc0-cc0_n)/cc0,abs(cc1-cc1_n)/cc1, abs(cc2-cc2_n)/cc2, abs(cc3-cc3_n)/cc3,...
    abs(va-va_n)/abs(va), abs(vc-vc_n)/abs(vc), abs(xx-xx_n)/xx];
err_cell_norm(n)=max(err_cell_vec); % maximum norm
% err_cell_norm(n)=vecnorm(err_cell_vec); % 2-norm

% relative cleft error: each row contains the error on the concentration of
% 1 solute/the potential/ pressure
err_cleft_mat=[abs(cl0-cl0_n)/vecnorm(cl0);abs(cl1-cl1_n)/vecnorm(cl1);abs(cl2-cl2_n)/vecnorm(cl2);...
    abs(cl3-cl3_n)/vecnorm(cl3);abs(vl-vl_n)/vecnorm(vl);abs(p_n-p)/vecnorm(p)]; 
 
row_norm=zeros(1,6);
for i=1:6
   % row_norm(i)=vecnorm(err_cleft_mat(i,:)); % 2-norm
   row_norm(i)=max(err_cleft_mat(i,:)); % max norm
end
err_cleft_norm(n)=max(row_norm);
% err_cleft_norm(n)=norm(err_cleft_mat,Inf);

% bvp5c error 
err_bvp(n)=bvperr_n;

%% update the variables

% cleft
cl0=cl0_n; cl1=cl1_n; cl2=cl2_n; cl3=cl3_n; vl=vl_n; dvl=dvl_n;
dcl0=dcl0_n; dcl1=dcl1_n; dcl2=dcl2_n; p=p_n; dp=dp_n;
y1=y1_n; y2=y2_n; y3=y3_n; y4=y4_n; y5=y5_n; y6=y6_n; y7=y7_n; y8=y8_n; y9=y9_n; y10=y10_n;

% cell
cc0=cc0_n;cc1=cc1_n; cc2=cc2_n; cc3=cc3_n; vc=vc_n; xx=xx_n; va=va_n;

% fluxes
J0bc=J0bc_n; J1bc=J1bc_n; J2bc=J2bc_n; J3bc=J3bc_n; 
J0ca=J0ca_n; J1ca=J1ca_n; J2ca=J2ca_n; J3ca=J3ca_n; 

%% check the residuals for the cell equations
% (exactly those solved, because after the cell the varibles were not
% updated again)

% non-dimensional potentials
Phi_ca= vc-va;
Phi_lc= vl-vc*ones(1,npts); % vector
Phi_la= vl(npts)-va; 

Jkira=C*Pkira*Phi_ca*(cc1 - ca1*exp(-Phi_ca))/(1 - exp(-Phi_ca));
Jcftr=-C*Pcftr*Phi_ca*(cc2 - ca2*exp(Phi_ca))/(1 - exp(Phi_ca));
Jn3b = Pn3b*(log(cc0*cc3^3/ca0/ca3^3)-2*Phi_ca);

% dissociation constants for Na-K pump
K0=0.2*(1+cc1*C/8.33); % dimensional
K0st=K0/C; % dimensionless

%%% pointwise ion fluxes lateral to cell (VECTORS) through the lateral membrane

K1_l= 0.1*(1 + cl0*C/(18.5));
K1_l_st=K1_l/C;

Jpump_l = Ppump*(cc0/(cc0 + K0st))^3*(cl1./(cl1 + K1_l_st)).^2;
Jn2b_l = Pn2b*(log(cl0.*cl3.^2/cc0/cc3^2)-Phi_lc);% + out of the cell
Jnkcc_l = Pnkcc*log(cl0.*cl1.*cl2.^2/cc0/cc1/cc2^2);
Jae_l = Pae*log(cl2.*cc3/cc2./cl3); % basolateral
Jkirb_l = Pkirb*C*Phi_lc.*(cl1 - cc1.*exp(-Phi_lc))./(1 - exp(-Phi_lc));

J0lc=-3*Jpump_l+Jnkcc_l+Jn2b_l;
J1lc=2*Jpump_l+Jnkcc_l+Jkirb_l;
J2lc=2*Jnkcc_l+Jae_l;
J3lc=2*Jn2b_l-Jae_l;

%%% ion flux across tj, considering as lateral values those at the lateral 
% side of the tj (x=L)
% include area factors here in the fluxes
Jtj=Atj*Ptj*(Phi_la)*C*(cl0(npts) - ca0*exp(-Phi_la))/(1-exp(-Phi_la));
J1tj=Atj*P1tj*(Phi_la)*C*(cl1(npts) - ca1*exp(-Phi_la))/(1-exp(-Phi_la));
J2tj=-Atj*P2tj*(Phi_la)*C*(cl2(npts) - ca2*exp(Phi_la))/(1-exp(Phi_la));
J3tj=-Atj*P3tj*(Phi_la)*C*(cl3(npts) - ca3*exp(Phi_la))/(1-exp(Phi_la));
    
% dimensionless equations
% ions
eq0=J0bc/Ppump+2*L*sum(w_cheb_01.*J0lc)/w/Ppump-J0ca/Ppump;
eq1=J1bc/Ppump+2*L*sum(w_cheb_01.*J1lc)/w/Ppump-J1ca/Ppump;
eq2=J2bc/Ppump+2*L*sum(w_cheb_01.*J2lc)/w/Ppump-J2ca/Ppump;
eq3=J3bc/Ppump+2*L*sum(w_cheb_01.*J3lc)/w/Ppump-J3ca/Ppump;

% water balance
uy_nd=xx+cc0+cc1+cc2+cc3-cl0-cl1-cl2-cl3;
eq_wb = (xx+cc0+cc1+cc2+cc3-cb0-cb1-cb2-cb3)+2*L*sum(w_cheb_01.*uy_nd)/w+...
    2*L^2*mu*Ux*sum(w_cheb_01.*p)/w/h^2/C/R/T-...
    (ca0+ca1+ca2+ca3-cc0-cc1-cc2-cc3-xx);

% electroneutrality cell
eq_en = cc0 +cc1 -cc2-cc3+zx*xx ; 

% open circuit
eq_oc= (Jtj/Ptj/C+J1tj/Ptj/C-J2tj/Ptj/C-J3tj/Ptj/C)+(-Jcftr/Ptj/C+Jkira/Ptj/C-2*Jn3b/Ptj/C);
% don't need area factors here, it's in the tj fluxes

% residuals
res_cell(n)=max(abs([eq0;eq1;eq2;eq3;eq_wb;eq_en;eq_oc]));


%% stopping criteria
if err_cleft_norm(n)<toll && err_cell_norm(n)<toll
break
end

% if res_cell(n)<toll && bvperr_n<toll
% break
% end

% if err_cleft_norm(n)<toll && res_cell(n)<toll
% break
% end
end

%%%%%%%%%%%%%% out of the itertion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% here we have some tests to check the solution in the cleft
%% check boundary conditions for the cleft

bar_sigma=chi*(y7-vc)-bar_sigma0;

% water velocity from lateral to cell, osmotic + hydraulic pressure
v_osm_lc=Lp*(R*T*C*L*(xx+cc0+cc1+cc2+cc3-2*(y1+y3))/Ux/h+mu*L^2*y9/h^3);

% \hat{C} and its derivative
CofX=2*(y1+y3);
dCofX=2*(y2+y4);
% slip velocity
u_slip=-bar_sigma.*y8./sqrt(CofX)-bar_sigma.^2.*dCofX./8./CofX.^2;
% second derivative of the potential in the cleft
d2V=(-(y2+y4).*y8-half_fac*v_osm_lc.*(y1/D0+y3/D1-y5/D2-(y1+y3-y5)/D3)+...
    half_fac*H*(-H^2*y10/3+u_slip).*(y2/D0+y4/D1-y6/D2-(y2+y4-y6)/D3)+L^2*(J0lc/D0+J1lc/D1-J2lc/D2-J3lc/D3)/h/C)./(y1+y3);
% second derivative of the sodium concentration in the cleft
d2c0=-y2.*y8-y1.*d2V+fac*(-y1.*v_osm_lc-H^3*y2.*y10/3+H*y2.*u_slip)/D0+2*L^2*J0lc/h/C/D0;
% second derivative of the potassium concentration in the cleft
d2c1=-y4.*y8-y3.*d2V+fac*(-y3.*v_osm_lc-H^3*y4.*y10/3+H*y4.*u_slip)/D1+2*L^2*J1lc/h/C/D1;
% second derivative of the bicarbonate concentration in the cleft
d2c2=+y6.*y8+y5.*d2V+fac*(-y5.*v_osm_lc-H^3*y6.*y10/3+H*y6.*u_slip)/D2+2*L^2*J2lc/h/C/D2;
% second derivative of \hat{C} 
d2CofX=2*(d2c0+d2c1);

% water velocity from lateral to apical (tj), osmotic + hydraulic pressure
v_osm_tj=Ltj*(R*T*C*(ca0+ca1+ca2+ca3-2*(y1(end)+y3(end)))/Ux+mu*L*y9(end)/h^2);

% tight juntions without area factors here, I account for the area of the
% tj when I integrate the 2D expression 

% residuals for the cleft boundary conditions
res_bc_cleft=[y1(1)-cb0;...
    y3(1)-cb1;...
    y5(1)-cb2;...
    -y2(end)-y1(end)*y8(end)+half_fac*y1(end)*(-H^2*y10(end)/3+u_slip(end))/D0-Jtj*L/D0/C/Atj;...
    -y4(end)-y3(end)*y8(end)+half_fac*y3(end)*(-H^2*y10(end)/3+u_slip(end))/D1-J1tj*L/D1/C/Atj;...
    -y6(end)+y5(end)*y8(end)+half_fac*y5(end)*(-H^2*y10(end)/3+u_slip(end))/D2-J2tj*L/D2/C/Atj;...
    -(y2(end)+y4(end)-y6(end))+(y1(end)+y3(end)-y5(end))*y8(end)+half_fac...
    *(y1(end)+y3(end)-y5(end))*(-H^2*y10(end)/3+u_slip(end))-J3tj*L/D3/C/Atj;
    y7(1);
    y9(1);
    y10(end)-3*(u_slip(end)-v_osm_tj)/H^2];

%% check cleft equations

% derivative of u_slip
duslip=-chi.*y8.^2./sqrt(CofX)+bar_sigma.*dCofX.*y8./2./CofX./sqrt(CofX)-bar_sigma.*d2V./sqrt(CofX)...
    -bar_sigma.*chi.*y8.*dCofX./4./CofX.^2+bar_sigma.^2.*dCofX.^2./4./CofX.^3-bar_sigma.^2.*d2CofX./8./CofX.^2;

% second derivative of pressure
d2P=v_osm_lc*3/H^3+3*duslip/H.^2;

% residuals of the cleft equations
res_cleft=[max(D*y1'-y2'),...
    max(D*y2'-d2c0'),...
    max(D*y3'-y4'),...
    max(D*y4'-d2c1'),...
    max(D*y5'-y6'),...
    max(D*y6'-d2c2'),...
    max(D*y7'-y8'),...
    max(D*y8'-d2V'),...
    max(D*y9'-y10'),...
    max(D*y10'-d2P')];


%% check water conservation in the cleft

% dimensional concentrations
Cb0=C*cb0; Cb1=C*cb1; Cb2=C*cb2; Cb3=C*cb3; 
Cc0=C*cc0; Cc1=C*cc1; Cc2=C*cc2; Cc3=C*cc3; 
Ca0=ca0*C; Ca1=ca1*C; Ca2=ca2*C; Ca3=ca3*C; Xx=xx*C;

% dimensional poential
Va=va*R*T/F; Vc=vc*R*T/F;

% dimensional cell to apical flux across the cell membrane (apical)
Qca=w*Lp*R*T*(Ca0+Ca1+Ca2+Ca3-Cc0-Cc1-Cc2-Cc3-Xx);
% water flux across the tight junction from lateral to apical
qtj=-2*H^3*dp(end)/3+2*u_slip(end)*H; % dimensionless
Qtj=qtj*Ux*h; % dimensional m^2/s
% water flux at the inlet of the cleft
qin=-2*H^3*dp(1)/3+2*u_slip(1)*H; % dimensionless
Qin=qin*Ux*h; % dimensional m^2/s
% water flux across the lateral cell membrane, from cleft to cell
qlat=sum(w_cheb_01.*v_osm_lc); % dimensionless
Qlat=qlat*Ux*h; % dimensional m^2/s % dimensional

overall_water_flux=(qtj*Ux*h+Qca)/(w+h); % dimensional transendothelial water
% flux towards the apical compartment, in m/s

% check conservation of water (dimensionless)
Water_cons=qin-2*qlat-qtj; 


%% check ion conservations (dimensionless)

J_na_in_nd=-(dcl0(1)+dvl(1)*cl0(1))+fac*H*cl0(1)*(-H^2*dp(1)/3+u_slip(1))/D0;
J_na_out_nd=-(dcl0(end)+dvl(end)*cl0(end))+fac*H*cl0(end)*(-H^2*dp(end)/3+u_slip(end))/D0;
J_na_lat_nd=2*L^2*sum(w_cheb_01.*J0lc)/h/D0/C;
% check sodium conservation in the cleft
% J_na_in_nd-J_na_lat_nd-J_na_out_nd

J_k_in_nd=-(dcl1(1)+dvl(1)*cl1(1))+fac*H*cl1(1)*(-H^2*dp(1)/3+u_slip(1))/D1;
J_k_out_nd=-(dcl1(end)+dvl(end)*cl1(end))+fac*H*cl1(end)*(-H^2*dp(end)/3+u_slip(end))/D1;
J_k_lat_nd=2*L^2*sum(w_cheb_01.*J1lc)/h/D1/C;
% check potassium conservation in the cleft
% J_k_in_nd-J_k_lat_nd-J_k_out_nd

J_cl_in_nd=-(dcl2(1)-dvl(1)*cl2(1))+fac*H*cl2(1)*(-H^2*dp(1)/3+u_slip(1))/D2;
J_cl_out_nd=-(dcl2(end)-dvl(end)*cl2(end))+fac*H*cl2(end)*(-H^2*dp(end)/3+u_slip(end))/D2;
J_cl_lat_nd=2*L^2*sum(w_cheb_01.*J2lc)/h/D2/C;
% check potassium conservation in the cleft
% J_cl_in_nd-J_cl_lat_nd-J_cl_out_nd

% integrated fluxes
J0lc_int=2*L*sum(w_cheb_01.*J0lc)/w;
J1lc_int=2*L*sum(w_cheb_01.*J1lc)/w;
J2lc_int=2*L*sum(w_cheb_01.*J2lc)/w;
J3lc_int=2*L*sum(w_cheb_01.*J3lc)/w;

%%%%%%%%%%%%%% plots %%%%%%%%%%%%%%

%% error plots

% figure
% plot(log10(res_cell))
% title('log10 residuals cell')

% figure
% plot(log10(err_cell_norm))
% title('Cell error (log10)')

% figure
% plot(log10(err_cleft_norm))
% title('Cleft error (log10)')

% figure
% plot(log10(err_bvp))
% title('bvp error (log10)')

%% potential
figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
plot(xmesh*L,vl*R*T/F,'LineWidth',3)
%title('Potential along the cleft')
xlabel('$x$ (m)','Interpreter','Latex')
% ylabel('{\it V^l} (V)','FontName','Cambria Math')
ylabel('$V^l$ (V)','Interpreter','Latex')
xlim([0,xmesh(end)*L])
set(gca,'linewidth',3,'FontSize',28)
% saveas(gca,'Potential.fig')
% print(gcf, 'V_cleft.png', '-dpng', '-r200');


%% concentrations

% osmolarity in the cleft
figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
plot(xmesh*L,(cl0+cl1+cl2+cl3)*C,'LineWidth',3,'DisplayName','Cleft')
% title('Cleft osmolarity')
legend('Location','southwest')
xlabel('$x$ (m)','Interpreter','latex')
ylabel('Osmolarity (mol/m$^3$)','Interpreter','latex')
xlim([0,xmesh(end)*L])
hold on
plot(xmesh*L,ones(length(xmesh),1)*(cc0+cc1+cc2+cc3+xx)*C,'LineWidth',3,'DisplayName','Cell')
set(gca,'linewidth',3, 'FontSize',28)
% saveas(gca,'CleftCell Osmolarity')
% print(gcf, 'osm_cleft_cell.png', '-dpng', '-r200');

%% cleft vs basal concentrations

figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
plot(xmesh*L,(cl0-cb0)*C,'LineWidth',3,'DisplayName','Na$^+$')
% title('Concentrations in the cleft')
legend('Location','southwest')
xlabel('$x$ (m)','Interpreter','latex')
ylabel('Concentration (mol/m$^3$)','Interpreter','latex')
xlim([0,xmesh(end)*L])
hold on
plot(xmesh*L,(cl1-cb1)*C,'LineWidth',3,'DisplayName','K$^+$')
plot(xmesh*L,(cl2-cb2)*C,'LineWidth',3,'DisplayName','Cl$^-$')
plot(xmesh*L,(cl3-cb3)*C,'LineWidth',3,'DisplayName','HCO$_3^-$')
% plot(xmesh(1)*L,(cb0-cb0)*C,'*','LineWidth',15,'DisplayName','Basal')
% plot(xmesh(end)*L,(cb0-cb0)*C,'*','LineWidth',15,'DisplayName','Apical')
set(gca,'linewidth',3, 'FontSize',28)
% saveas(gca,'Cleft vs Basal concentrations')
% print(gcf, 'concVsBasal.png', '-dpng', '-r200');


%% slip velocity 

figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
plot(xmesh*L,u_slip*Ux,'LineWidth',3)
xlabel('$x$ (m)','Interpreter','latex')
xlim([xmesh(1) xmesh(end)]*L)
ylabel('Slip velocity (m/s)','Interpreter','latex')
set(gca,'linewidth',3, 'FontSize',28)
% saveas(gca,'slip_velocity')
% print(gcf, 'slip_velocity.png', '-dpng', '-r200');

%% velocity plots
%% ux

[w_cheb,y_cheb]=clenshaw_curtis_p(npts,1); % y_cheb in [1,-1]
ymesh=flip(x_cheb)'; % flip [-1,1]
ymesh_H=ymesh*H;
w_cheb=flip(w_cheb'); % weights for integration formula in -1 1
w_cheb_H=w_cheb*H; % weights for integration formula in -H H
ux=zeros(npts);
Qx=zeros(1,npts);
for i=1:npts
    y=ymesh_H(i);
    ux(i,:)=0.5*y10*(y^2-H^2)+u_slip;
end

for i=1:npts
    Qx(i)=sum(w_cheb_H'.*ux(:,i));
end

% figure
% surf(xmesh*L,ymesh_H*h,ux*Ux)
% xlabel('x (m)')
% ylabel('y (m)')
% title('u_x (m/s)')
% set(gca,'linewidth',2, 'FontSize',24, 'FontWeight','bold')

%% uy

uy=zeros(npts);
for i=1:npts
    y=ymesh_H(i);
    uy(i,:)=-0.5*d2P*y*(y^2/3-H^2)-duslip*y;
end

% figure
% surf(xmesh*L,ymesh_H*h,uy*Ux*h/L)
% xlabel('x (m)')
% ylabel('y (m)')
% % title('u_y (m/s)')
% set(gca,'linewidth',2, 'FontSize',24, 'FontWeight','bold')

%% slip velocity
figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
yyaxis left
plot(xmesh*L,Qx*Q,'LineWidth',3)
xlabel('$x$ (m)','Interpreter','latex')
xlim([xmesh(1) xmesh(end)]*L)
ylabel('Water flux (m$^2$/s)','Interpreter','latex')
yyaxis right
plot(xmesh*L,u_slip*Ux,'LineWidth',3)
ylabel('Slip velocity (m/s)','Interpreter','latex')
set(gca,'linewidth',3, 'FontSize',28)
% saveas(gca,'Water flux slip')
% print(gcf, 'water_flux_slip.png', '-dpng', '-r200');

%% water flux from the cleft to the cell

figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
plot(xmesh*L,Qx*Q,'LineWidth',3)
xlabel('$x$ (m)','Interpreter','latex')
xlim([xmesh(1) xmesh(end)]*L)
ylabel('Water flux (m$^2$/s)','Interpreter','latex')
set(gca,'linewidth',3, 'FontSize',28)
% saveas(gca,'Water flux')
% print(gcf, 'water_flux.png', '-dpng', '-r200');

%% contourplot 

figure('WindowState','maximized')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
mod_u=sqrt((ux*Ux).^2+(uy*Ux*h/L).^2);
cont_plot=pcolor(xmesh,ymesh_H,mod_u);
cont_plot.LineWidth=3;
col=colorbar;
col.Label.Interpreter = 'latex';
col.Label.String = '$|{\bf u}|$, m/s';
col.TickLabelInterpreter="latex";
shading interp
hold on
quiver(xmesh(1:10:end),ymesh_H(1:10:end),ux(1:10:end,1:10:end)*Ux,...
    uy(1:10:end,1:10:end)*Ux,'AutoScale','on','AutoScaleFactor',0.5,'LineWidth',2,'Color','k')
xlabel('$x/L$ ','Interpreter','latex')
ylabel('$y/h$ ','Interpreter','latex')
xlim([xmesh(1) xmesh(end)]);
ylim([ymesh_H(1) ymesh_H(end)]);

set(gca,'linewidth',2, 'FontSize',28)
% saveas(gca,'velocity')
% print(gcf, 'velocity.png', '-dpng', '-r200');

%% ion fluxes across the lateral membrane
% figure
% plot(xmesh*L,J0lc,'r','LineWidth',2,'DisplayName','Na')
%  hold on
% plot(xmesh*L,J1lc,'b','LineWidth',2,'DisplayName','K')
% plot(xmesh*L,J2lc,'g','LineWidth',2,'DisplayName','Cl')
% plot(xmesh*L,J3lc,'m','LineWidth',2,'DisplayName','HCO3')
%  hold off
% title('Ion fluxes across lateral membrane (lateral to cell)')
% legend
% xlabel('m')
% ylabel('m/s')
% xlim([0,xmesh(end)*L])
% set(gca,'linewidth',1.5, 'FontSize',12, 'FontWeight','bold')
% saveas(gca,'Ion_fluxes.fig')

%% other concentration plots
% figure
% plot(xmesh*L,cl0*C,'r','LineWidth',2,'DisplayName','Na cleft')
% %title('Na concentration in the cleft. Na in the cell is',cc0*C)
% legend
% xlabel('m')
% ylabel('mol/m^3')
% xlim([0,xmesh(end)*L])
% set(gca,'linewidth',2, 'FontSize',28, 'FontWeight','bold')
% saveas(gca,'Na_cleft')

% figure
% plot(xmesh*L,cl1*C,'b','LineWidth',3,'DisplayName','K')
% %title('K concentration in the cleft. K in the cell is',cc1*C)
% legend
% xlabel('m')
% ylabel('mol/m^3')
% xlim([0,xmesh(end)*L])
% set(gca,'linewidth',3, 'FontSize',28, 'FontWeight','bold')
% % saveas(gca,'K_cleft')

% figure
% plot(xmesh*L,cl2*C,'g','LineWidth',3,'DisplayName','Cl')
% %title('Cl concentration in the cleft. Cl in the cell is',cc2*C)
% legend
% xlabel('m')
% ylabel('mol/m^3')
% xlim([0,xmesh(end)*L])
% set(gca,'linewidth',3, 'FontSize',28, 'FontWeight','bold')
% saveas(gca,'Cl_cleft')

% figure
% plot(xmesh*L,cl3*C,'m','LineWidth',3,'DisplayName','HCO3')
% %title('HCO3 concentration in the cleft. HCO3 in the cell is',cc3*C)
% legend
% xlabel('m')
% ylabel('mol/m^3')
% xlim([0,xmesh(end)*L])
% set(gca,'linewidth',3, 'FontSize',28, 'FontWeight','bold')
% saveas(gca,'HCO3_cleft')

% figure
% plot(xmesh*L,cl1*C-Cb1,'g','LineWidth',2,'DisplayName','K^+')
% hold on
% plot(xmesh*L,cl2*C-Cb2,'b','LineWidth',2,'DisplayName','Cl^-')
% plot(xmesh*L,cl3*C-Cb3,'c','LineWidth',2,'DisplayName','HCO_3^-')
% plot(xmesh*L,cl0*C-Cb0,'k','LineWidth',2,'DisplayName','Na^+')
% title('Cleft concentrations with respect to the basal compartment')
% legend
% xlabel('m')
% ylabel('mol/m^3')
% xlim([0,xmesh(end)*L])
% set(gca,'linewidth',1.5, 'FontSize',12, 'FontWeight','bold')
% 
% figure
% plot(xmesh*L,(cl0+cl1+cl2+cl3)*C-(Cc0+Cc1+Cc2+Cc3+Xx),'k','LineWidth',2,'DisplayName','Osmolarity')
% title('Osmolarity in the cleft with respect to the cell')
% xlabel('m')
% ylabel('mol/m^3')
% xlim([0,xmesh(end)*L])
% set(gca,'linewidth',1.5, 'FontSize',24, 'FontWeight','bold')
% 
% figure(15)
% plot(xmesh*L,(cc0+cc1+cc2+cc3+xx)*C*ones(length(xmesh),1)-(cb0+cb1+cb2+cb3)*C, 'c','LineWidth',2,'DisplayName','Cell' )
% hold on
% plot(xmesh*L,(cl0+cl1+cl2+cl3)*C-(cb0+cb1+cb2+cb3)*C,'k','LineWidth',2,'DisplayName','Cleft')
% title('Osmolarity in the cleft and in the cell with respect to basal')
% legend
% % % xlabel('m')
% ylabel('mol/m^3')
% xlim([0,xmesh(end)*L])
% set(gca,'linewidth',1.5, 'FontSize',12, 'FontWeight','bold')
% 

%% pressure
osm_press_lc_dim= R*T*(cc0+cc1+cc2+cc3+xx-cl0-cl1-cl2-cl3)*C;
osm_press_la_dim= R*T*(ca0+ca1+ca2+ca3-cl0(end)-cl1(end)-cl2(end)-cl3(end))*C;

mech_press_lc_dim=y9*mu*Ux*L/h^2;
% figure
% plot(xmesh*L,y9*mu*Ux*L/h^2,'LineWidth',3)
% title('Dimensional pressure')
% xlim([0,xmesh(end)*L])
% xlabel('m')
% ylabel('Pa')
% set(gca,'linewidth',3, 'FontSize',28, 'FontWeight','bold')
% 
% hold on
% plot(xmesh*L,osm_press_lc_dim,'LineWidth',3)
% plot(xmesh(end)*L,osm_press_la_dim,'-o')
% legend('mechanical pressure', 'osmotic pressure cleft-cell','osmotic pressure at the tj','Location','northwest')
% ylabel('Pa')
% saveas(gca,'Pressure')


%% area factors to compute water flux across endothelium  
% hexagonal cells 

Acell=(w/2)^2*3*sqrt(3)/2;
Atot=(w/2+h/sqrt(3))^2*3*sqrt(3)/2;
Acleft=Atot-Acell;

uTJ=Qtj/h; % velocity of water at the tight junction
uca=Qca/w; % velocity of water at the apical cell membrane
% mean velocity across the endothelium accounting for the area factors
Qcell=uca*Acell;
Qcleft=uTJ*Acleft;
u_end_hex=1/Atot*(Qcell+Qcleft);

%% display
% cellular concentrations
fprintf('Sodium concentration in the cell = %g (mM)\n',Cc0)
fprintf('Potassium concentration in the cell = %g (mM)\n',Cc1)
fprintf('Chloride concentration in the cell = %g (mM)\n',Cc2)
fprintf('Bicarbonate concentration in the cell = %g (mM)\n',Cc3)
% electrical potential
fprintf('Transendothelial potential = %g (V)\n',Va)
fprintf('Cellular potential = %g (V)\n',Vc)
% water fluxes
fprintf('Transendothelial water flux = %g m/s\n',overall_water_flux)
fprintf('Cell-to-apical water flux = %g m^2/s\n',Qca)
fprintf('Water flux across the tj = %g m^2/s\n',Qtj)
fprintf('Water flux across the lateral membrane = %g m^2/s\n',2*Qlat)
fprintf('Water flux at the inlet of the cleft = %g m^2/s\n',Qin)
% transendothelial ion fluxes
fprintf('Transendothelial sodium flux = %g (mol/m^2/s)\n',J0ca+Jtj)
fprintf('Transendothelial potassium flux = %g (mol/m^2/s)\n',J1ca+J1tj)
fprintf('Transendothelial chloride flux = %g (mol/m^2/s)\n',J2ca+J2tj)
fprintf('Transendothelial bicarbonate flux = %g (mol/m^2/s)\n',J3ca+J3tj)
% tight junction ion fluxes
fprintf('Tight junction sodium flux = %g (mol/m^2/s)\n',Jtj)
fprintf('Tight junction potassium flux = %g (mol/m^2/s)\n',J1tj)
fprintf('Tight junction chloride flux = %g (mol/m^2/s)\n',J2tj)
fprintf('Tight junction bicarbonate flux = %g (mol/m^2/s)\n',J3tj)
% integrated flux across the lateral membrane
fprintf('Lateral sodium flux = %g (mol/m^2/s)\n',J0lc_int)
fprintf('Lateral potassium flux = %g (mol/m^2/s)\n',J1lc_int)
fprintf('Lateral chloride flux = %g (mol/m^2/s)\n',J2lc_int)
fprintf('Lateral bicarbonate flux = %g (mol/m^2/s)\n',J3lc_int)









