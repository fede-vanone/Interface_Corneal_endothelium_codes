function obj= objective_nd(permeabilities)
% this function computes the objective for the optimization in the cellular
% domain, given a vector of 13 permeabilities as input

Ppump_nd=permeabilities(1);
Pnkcc_nd=permeabilities(2);
Pkirb_nd=permeabilities(3);
Pkira_nd=permeabilities(4);
Pcftr_nd=permeabilities(5);
Pn2b_nd=permeabilities(6);
Pn3b_nd=permeabilities(7);
Paea_nd=permeabilities(8);
Pae_nd=permeabilities(9);
P0tj_nd=permeabilities(10);
P1tj_nd=permeabilities(11);
P2tj_nd=permeabilities(12);
P3tj_nd=permeabilities(13);

% dimensionless parameters
[F,R,T,h,L,~,~,C,zx,Lp,Ltj,...
    Atj,w,~,~,~,~,Cb0_nd,Cb1_nd,Cb2_nd,Cb3_nd,Vb_nd,Ca0_nd,Ca1_nd,Ca2_nd,Ca3_nd,...
    toll,nit,npts]=pars5();

% define symbolic variables
syms Cc0_nd Cc1_nd Cc2_nd Cc3_nd Xx_nd Vc_nd Va_nd 
% 0 Na
% 1 K
% 2 Cl
% 3 HCO3

%%
% dissociation constants, dimensionless
Kna_nd = 0.2*(1 + Cc1_nd*C/(8.33))/C; Kk_nd = 0.1*(1 + Cb0_nd*C/(18.5))/C;

% dimensionless potentials
Phi_bc= Vb_nd-Vc_nd;
Phi_ca= Vc_nd-Va_nd;
Phi_ba= Vb_nd-Va_nd;

%% fluxes

% active
Jpump_nd = Ppump_nd*(Cc0_nd/(Cc0_nd + Kna_nd))^3*(Cb1_nd/(Cb1_nd + Kk_nd))^2;
Jnkcc_nd = Pnkcc_nd*log(Cb0_nd*Cb1_nd*Cb2_nd^2/Cc0_nd/Cc1_nd/Cc2_nd^2);

Jn2b_nd = Pn2b_nd*(log(Cb0_nd*Cb3_nd^2/Cc0_nd/Cc3_nd^2)-Phi_bc);% + out of the cell
Jn3b_nd = Pn3b_nd*(log(Cc0_nd*Cc3_nd^3/Ca0_nd/Ca3_nd^3)-2*Phi_ca);% + out of the cell

Jae_nd = Pae_nd*log(Cb2_nd*Cc3_nd/Cc2_nd/Cb3_nd); % basolateral
Jaea_nd=Paea_nd*log(Ca2_nd*Cc3_nd/Cc2_nd/Ca3_nd); % apical

% passive
Jkirb_nd = Pkirb_nd*Phi_bc*(Cb1_nd - Cc1_nd*exp(-Phi_bc))/(1 - exp(-Phi_bc));
Jkira_nd=Pkira_nd*Phi_ca*(Cc1_nd - Ca1_nd*exp(-Phi_ca))/(1 - exp(-Phi_ca));
Jcftr_nd=-Pcftr_nd*Phi_ca*(Cc2_nd - Ca2_nd*exp(Phi_ca))/(1 - exp(Phi_ca));

% tj 
J0tj_nd=Atj*P0tj_nd*Phi_ba*(Cb0_nd - Ca0_nd*exp(-Phi_ba))/(1-exp(-Phi_ba));
J1tj_nd=Atj*P1tj_nd*Phi_ba*(Cb1_nd - Ca1_nd*exp(-Phi_ba))/(1-exp(-Phi_ba));
J2tj_nd=-Atj*P2tj_nd*Phi_ba*(Cb2_nd - Ca2_nd*exp(Phi_ba))/(1-exp(Phi_ba));
J3tj_nd=-Atj*P3tj_nd*Phi_ba*(Cb3_nd - Ca3_nd*exp(Phi_ba))/(1-exp(Phi_ba));

%% cell equations
% balance in the cell

eq0 = (1+2*L/w)*(-3*Jpump_nd + Jnkcc_nd + Jn2b_nd) -Jn3b_nd == 0; % Na
eq1 = (1+2*L/w)*(2*Jpump_nd + Jnkcc_nd +Jkirb_nd) - Jkira_nd == 0; % K
eq2 = (1+2*L/w)*(2*Jnkcc_nd + Jae_nd) - Jcftr_nd +Jaea_nd== 0; % Cl
eq3 = (1+2*L/w)*(2*Jn2b_nd -Jae_nd) -3*Jn3b_nd -Jaea_nd== 0; % HCO3

% water balance
eq_wb = (Ca0_nd+Ca1_nd+Ca2_nd+Ca3_nd-Cc0_nd-Cc1_nd-Cc2_nd-Cc3_nd-Xx_nd)-...
    (1+2*L/w)*(Xx_nd+Cc0_nd+Cc1_nd+Cc2_nd+Cc3_nd-Cb0_nd- Cb1_nd-Cb2_nd-Cb3_nd)==0;

% electroneutrality cell
eq_en = Cc0_nd +Cc1_nd -Cc2_nd-Cc3_nd+zx*Xx_nd ==0; 

% open circuit
eq_oc= J0tj_nd + J1tj_nd-J2tj_nd-J3tj_nd+(-Jcftr_nd+Jkira_nd-2*Jn3b_nd)==0; % area factor already inside Jtj


%% cell solution

tic
eqns=[eq0 eq1 eq2 eq3 eq_wb eq_en eq_oc];
vars = [Cc0_nd Cc1_nd Cc2_nd Cc3_nd Xx_nd Vc_nd Va_nd];
range = [0 100/C;0 200/C;0 100/C;0 100/C;0 200/C; -3.743601232278240 0;-0.374360123227824 0.374360123227824];
clear sol

sol = vpasolve(eqns,vars,range);
toc

%% output processing

% concentrations dimensionless
Cc0_nd=double(sol.Cc0_nd);
Cc1_nd=double(sol.Cc1_nd);
Cc2_nd=double(sol.Cc2_nd);
Cc3_nd=double(sol.Cc3_nd);

% potentials
Vc_nd=double(sol.Vc_nd);
Va_nd=double(sol.Va_nd);

%% objective 

Vc=Vc_nd*R*T/F;
Va=Va_nd*R*T/F;
% reference values 
Cc0_exp=15/C; 
Cc1_exp=132/C;
Cc2_exp=38/C;
Cc3_exp=25/C; 
Vc_exp=(-60e-3);
Va_exp=(-.5e-3);


obj=sqrt((Cc0_nd-Cc0_exp)^2/Cc0_exp^2+(Cc1_nd-Cc1_exp)^2/Cc1_exp^2+(Cc2_nd-Cc2_exp)^2/Cc2_exp^2+...
   (Cc3_nd-Cc3_exp)^2/Cc3_exp^2+(Vc-Vc_exp)^2/Vc_exp^2+(Va-Va_exp)^2/Va_exp^2);


end