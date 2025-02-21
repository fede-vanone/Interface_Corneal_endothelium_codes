% This code is for optimizing the channel permeabilities on the cellular domain, 
% targeting cellular concentrations cellular and transendothelial potential.
% 
% Once computed the permeabilities, the model is run (only on the cellular domain)
% with those as parameters.
% The optimized permeabilities are in 
% permeabilities_vec=[Ppump,Pnkcc,Pkirb,Pkira,Pcftr,Pn2b,Pn3b,Paea,Pae,P0tj,P1tj,P2tj,P3tj]

clear all
%% initialize the solver
[F,R,T,h,L,~,~,C,zx,Lp,Ltj,...
    Atj,w,~,~,~,~,Cb0_nd,Cb1_nd,Cb2_nd,Cb3_nd,Vb_nd,Ca0_nd,Ca1_nd,Ca2_nd,Ca3_nd,...
    toll,nit,npts]=pars5();

ref=1e-6;
% define the starting point for fminsearch
par_guess = 5*[10^(-6),5*10^(-7),.5*10^(-7)*C,.5*10^(-7)*C,3*10^(-8)*C,10^(-6),.5*10^(-6),...
    10^(-7),10^(-7),(6e-7)*C*w/h,(7.7e-7)*C*w/h,(3.2e-7)*C*w/h,(3.2e-7)*C*w/h]/ref;

%% solve with fminsearch 

par_vec=par_opt_nd(par_guess);
% display the value of the objective function at each iteration
objval=objective_nd(par_vec)

%% get permeabilities (dimensionless)

Ppump_nd=par_vec(1);
Pnkcc_nd=par_vec(2);
Pkirb_nd=par_vec(3);
Pkira_nd=par_vec(4);
Pcftr_nd=par_vec(5);
Pn2b_nd=par_vec(6);
Pn3b_nd=par_vec(7);
Paea_nd=par_vec(8);
Pae_nd=par_vec(9);
P0tj_nd=par_vec(10);
P1tj_nd=par_vec(11);
P2tj_nd=par_vec(12);
P3tj_nd=par_vec(13);


%% run the model with the optimized permeabilities

syms Cc0_nd Cc1_nd Cc2_nd Cc3_nd Xx_nd Vc_nd Va_nd 
% 0 Na
% 1 K
% 2 Cl
% 3 HCO3

%%
% dissociation constants, dimensionless
Kna_nd = 0.2*(1 + Cc1_nd*C/(8.33))/C; Kk_nd = 0.1*(1 + Cb0_nd*C/(18.5))/C;

% non-dimensional potentials
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


%% cellular equations 
% balance in the cell

eq0 = (1+2*L/w)*(-3*Jpump_nd + Jnkcc_nd + Jn2b_nd) -Jn3b_nd == 0; % Na
eq1 = (1+2*L/w)*(2*Jpump_nd + Jnkcc_nd +Jkirb_nd) - Jkira_nd == 0; % K
eq2 = (1+2*L/w)*(2*Jnkcc_nd + Jae_nd) - Jcftr_nd +Jaea_nd== 0; % Cl
eq3 = (1+2*L/w)*(2*Jn2b_nd -Jae_nd) -3*Jn3b_nd -Jaea_nd== 0; % HCO3


% water balance
eq_wb = (Ca0_nd+Ca1_nd+Ca2_nd+Ca3_nd-Cc0_nd-Cc1_nd-Cc2_nd-Cc3_nd-Xx_nd)-...
    (1+2*L/w)*(Xx_nd+Cc0_nd+Cc1_nd+Cc2_nd+Cc3_nd-Cb0_nd- Cb1_nd-Cb2_nd-Cb3_nd)==0;

% electroneutrality cell
eq_en = Cc0_nd + Cc1_nd -Cc2_nd-Cc3_nd+zx*Xx_nd ==0; 

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
Xx_nd=double(sol.Xx_nd);

% potentials
Vc_nd=double(sol.Vc_nd);
Va_nd=double(sol.Va_nd);


%% recompute dimensionless fluxes

% dissociation constants, dimensionless
Kna_nd = 0.2*(1 + Cc1_nd*C/(8.33))/C; Kk_nd = 0.1*(1 + Cb0_nd*C/(18.5))/C;

% non-dimensional potentials
Phi_bc= Vb_nd-Vc_nd;
Phi_ca= Vc_nd-Va_nd;
Phi_ba= Vb_nd-Va_nd;

%% check objective function
Cc0_exp=15/C; 
Cc1_exp=132/C;
Cc2_exp=38/C;
Cc3_exp=25/C; 
Vc_exp=(-60e-3)*F/R/T;
Va_exp=(-.5e-3)*F/R/T;

obj=sqrt((Cc0_nd-Cc0_exp)^2/Cc0_exp^2+(Cc1_nd-Cc1_exp)^2/Cc1_exp^2+(Cc2_nd-Cc2_exp)^2/Cc2_exp^2+...
   (Cc3_nd-Cc3_exp)^2/Cc3_exp^2+(Vc_nd-Vc_exp)^2/Vc_exp^2+(Va_nd-Va_exp)^2/Va_exp^2);

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

% dimensionalize fluxes 
Jpump=Jpump_nd*ref;
Jnkcc=Jnkcc_nd*ref;
Jn2b=Jn2b_nd*ref;
Jn3b=Jn3b_nd*ref;
Jae=Jae_nd*ref;
Jaea=Jaea_nd*ref;
Jkirb=Jkirb_nd*ref;
Jkira=Jkira_nd*ref;
Jcftr=Jcftr_nd*ref;
J0tj=J0tj_nd*ref;
J1tj=J1tj_nd*ref;
J2tj=J2tj_nd*ref;
J3tj=J3tj_nd*ref;

%% permeabilities

% dimensional
Ppump=Ppump_nd*ref;
Pnkcc=Pnkcc_nd*ref;
Pkirb=Pkirb_nd*ref/C;
Pkira=Pkira_nd*ref/C;
Pcftr=Pcftr_nd*ref/C;
Pn2b=Pn2b_nd*ref;
Pn3b=Pn3b_nd*ref;
Paea=Paea_nd*ref;
Pae=Pae_nd*ref;
P0tj=P0tj_nd*ref/C;
P1tj=P1tj_nd*ref/C;
P2tj=P2tj_nd*ref/C;
P3tj=P3tj_nd*ref/C;

permeabilities_vec=[Ppump,Pnkcc,Pkirb,Pkira,Pcftr,Pn2b,Pn3b,Paea,Pae,P0tj,P1tj,P2tj,P3tj]'

%% dimensional variables
Cc0=Cc0_nd*C; 
Cc1=Cc1_nd*C; 
Cc2=Cc2_nd*C; 
Cc3=Cc3_nd*C;
Vc=R*T*Vc_nd/F; 
Va=R*T*Va_nd/F;

% check the fluxes
J0ca=Jn3b;
J0bc=-3*Jpump + Jnkcc + Jn2b;
J0lateral=2*L*J0bc/w;
J0tj+J0ca; % transendothelial sodium flux

J1ca=Jkira;
J1bc=2*Jpump + Jnkcc +Jkirb;
J1lateral=2*L*J1bc/w;
J1tj+J1ca; % transendothelial potassium flux

J2ca=Jcftr-Jaea;
J2bc=2*Jnkcc + Jae;
J2lateral=2*L*J2bc/w;
J2tj+J2ca; % transendothelial chloride flux

J3ca=3*Jn3b +Jaea;
J3bc=2*Jn2b -Jae;
J3lateral=2*L*J3bc/w;
J3tj+J3ca; % transendothelial bicarbonate flux

