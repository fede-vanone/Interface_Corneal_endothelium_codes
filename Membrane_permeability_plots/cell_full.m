function [cc0,cc1,cc2,cc3,vc,va,xx,J0bc,J1bc,J2bc,J3bc,J0ca,J1ca,J2ca,J3ca]=cell_full(npts,cl0,cl1,cl2,cl3,vl,p)

global F R T Lp Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx Ptj Atj P1tj P2tj P3tj w...
 ca0 ca1 ca2 ca3 ...
 cb0 cb1 cb2 cb3 vb ...
 xmesh L h w_cheb_01 C mu Ux


%% initialization

% symbolic variables: cellular concentrations, cellular and apical
% potential
syms cc0 cc1 cc2 cc3 vc va xx 
% 0 Na
% 1 K
% 2 Cl
% 3 HCO3

%% fluxes

% dissociation constants for Na-K pump
K0=0.2*(1+cc1*C/8.33); % dimensional
K0st=K0/C; % dimensionless

K1=0.1*(1+cb0*C/18.5); % dimensional
K1st=K1/C; % dimensionless

% Kna = 0.2*(1 + Cc1/(8.38)); Kk = 0.1*(1 + Cb0/(18.5));

% non-dimensional potentials
Phi_bc= vb-vc;
Phi_ca= vc-va;
Phi_lc= vl-vc*ones(1,npts); % vector
Phi_la= vl(npts)-va; 

%%% ion flux across each channel/transporter on basal and apical membrane

% Na/K pump
Jpump = Ppump*(cc0/(cc0 + K0st))^3*(cb1/(cb1 + K1st))^2;
% Na/2bic cotransporter
Jn2b = Pn2b*(log(cb0*cb3^2/cc0/cc3^2)-Phi_bc); % + out of the cell
% Na/K/2Cl cotransporter
Jnkcc = Pnkcc*log(cb0*cb1*cb2^2/cc0/cc1/cc2^2);
% basolateral anion exchanger
Jae = Pae*log(cb2*cc3/cc2/cb3); 
% basolateral K channel
Jkirb = C*Pkirb*Phi_bc*(cb1 - cc1*exp(-Phi_bc))/(1 - exp(-Phi_bc));
% Na/3bic cotransporter
Jn3b = Pn3b*(log(cc0*cc3^3/ca0/ca3^3)-2*Phi_ca);% + out of the cell
% apical anion exchanger
Jaea=Paea*log(ca2*cc3/cc2/ca3);
% apical K channel
Jkira=C*Pkira*Phi_ca*(cc1 - ca1*exp(-Phi_ca))/(1 - exp(-Phi_ca));
% Cl channel
Jcftr=-C*Pcftr*(Phi_ca)*(cc2 - ca2*exp(Phi_ca))/(1 - exp(Phi_ca));

% overall ion fluxes basal to cell through the basal membrane
J0bc=-3*Jpump+Jnkcc+Jn2b;
J1bc=2*Jpump+Jnkcc+Jkirb;
J2bc=2*Jnkcc+Jae;
J3bc=2*Jn2b-Jae;

% overall ion fluxes cell to apical through the apical membrane
J0ca=Jn3b;
J1ca=Jkira;
J2ca=Jcftr-Jaea;
J3ca=3*Jn3b+Jaea;

%%% pointwise ion fluxes lateral to cell (VECTORS) through the lateral membrane

K1_l= 0.1*(1 + cl0*C/(18.5));
K1_l_st=K1_l/C;

% 
Jpump_l = Ppump*(cc0/(cc0 + K0st))^3*(cl1./(cl1 + K1_l_st)).^2;
Jn2b_l = Pn2b*(log(cl0.*cl3.^2/cc0/cc3^2)-Phi_lc);% + out of the cell
Jnkcc_l = Pnkcc*log(cl0.*cl1.*cl2.^2/cc0/cc1/cc2^2);
Jae_l = Pae*log(cl2*cc3/cc2./cl3); % basolateral
Jkirb_l = Pkirb*C*(Phi_lc).*(cl1 - cc1.*exp(-Phi_lc))./(1 - exp(-Phi_lc));

J0lc=-3*Jpump_l+Jnkcc_l+Jn2b_l;
J1lc=2*Jpump_l+Jnkcc_l+Jkirb_l;
J2lc=2*Jnkcc_l+Jae_l;
J3lc=2*Jn2b_l-Jae_l;

%%% ion flux across tj, considering as lateral values those at the lateral 
% side of the tj (x=L)
% 
% Jtj=Ptj*Phi_la(npts)*(Cl0(npts) - Ca0*exp(-Phi_la(npts)))/(1-exp(-Phi_la(npts)));
% J1tj=P1tj*Phi_la(npts)*(Cl1(npts) - Ca1*exp(-Phi_la(npts)))/(1-exp(-Phi_la(npts)));
% J2tj=-P2tj*Phi_la(npts)*(Cl2(npts) - Ca2*exp(Phi_la(npts)))/(1-exp(Phi_la(npts)));
% J3tj=-P3tj*Phi_la(npts)*(Cl3(npts) - Ca3*exp(Phi_la(npts)))/(1-exp(Phi_la(npts)));

Jtj=Atj*Ptj*(vl(end)-va)*C*(cl0(npts) - ca0*exp(-vl(end)+va))/(1-exp(-vl(end)+va));
J1tj=Atj*P1tj*(vl(end)-va)*C*(cl1(npts) - ca1*exp(-vl(end)+va))/(1-exp(-vl(end)+va));
J2tj=-Atj*P2tj*(vl(end)-va)*C*(cl2(npts) - ca2*exp(vl(end)-va))/(1-exp(vl(end)-va));
J3tj=-Atj*P3tj*(vl(end)-va)*C*(cl3(npts) - ca3*exp(vl(end)-va))/(1-exp(vl(end)-va));
% these fluxes include area factors.
    
%% dimensionless equations
% ions
eq0=J0bc/Ppump+2*L*sum(w_cheb_01.*J0lc)/w/Ppump-J0ca/Ppump==0;
eq1=J1bc/Ppump+2*L*sum(w_cheb_01.*J1lc)/w/Ppump-J1ca/Ppump==0;
eq2=J2bc/Ppump+2*L*sum(w_cheb_01.*J2lc)/w/Ppump-J2ca/Ppump==0;
eq3=J3bc/Ppump+2*L*sum(w_cheb_01.*J3lc)/w/Ppump-J3ca/Ppump==0;

% water balance
% uy=Lp*R*T*C*(xx+cc0+cc1+cc2+cc3-cl0-cl1-cl2-cl3);
uy_nd=xx+cc0+cc1+cc2+cc3-cl0-cl1-cl2-cl3;
eq_wb = (xx+cc0+cc1+cc2+cc3-cb0-cb1-cb2-cb3)+2*L*sum(w_cheb_01.*uy_nd)/w+...
    2*L^2*mu*Ux*sum(w_cheb_01.*p)/w/h^2/C/R/T-...
    (ca0+ca1+ca2+ca3-cc0-cc1-cc2-cc3-xx)==0;

% electroneutrality cell
eq_en = cc0 +cc1 -cc2-cc3+zx*xx ==0; 

% open circuit
% eq_oc= h*(Jtj+J1tj-J2tj-J3tj)+w*(-Jcftr+Jkira-2*Jn3b)==0;
eq_oc= (Jtj/Ptj/C+J1tj/Ptj/C-J2tj/Ptj/C-J3tj/Ptj/C)+(-Jcftr/Ptj/C+Jkira/Ptj/C-2*Jn3b/Ptj/C)==0;
% don't need area factor here, it's in the fluxes


%% solution
eqns=[eq0 eq1 eq2 eq3 eq_wb eq_en eq_oc];
vars = [cc0 cc1 cc2 cc3 xx vc va];
% incon = [15/C,130/C,40/C,30/C,60/C,-65/1000*F/R/T,-0.5/1000*F/R/T];
% range = [0 10;0 20;0 10;0 10;0 20; -10 0;-10 0];
% range = [0 100/C;0 200/C;0 100/C;0 100/C;0 200/C; -10*F/R/T 10*F/R/T;-10*F/R/T 10*F/R/T];
% range = [0 100;0 200;0 100;0 100;0 200; -10 10;-10 10];
range = [0 100/C;0 200/C;0 100/C;0 100/C;0 200/C; -.1*F/R/T 0;-0.01*F/R/T 0.01*F/R/T];
clear sol

sol = vpasolve(eqns,vars,range);


%% output
% concentrations
cc0=double(sol.cc0);
cc1=double(sol.cc1);
cc2=double(sol.cc2);
cc3=double(sol.cc3);
xx=double(sol.xx);

% potentials
vc=double(sol.vc);
va=double(sol.va);
% TEP=Vs-Vp

%% 
% dissociation constants for Na-K pump
K0=0.2*(1+cc1*C/8.33); % dimensional
K0St=K0/C; % dimensionless

K1=0.1*(1+cb0*C/18.5); % dimensional
K1St=K1/C; % dimensionless

% Kna = 0.2*(1 + Cc1/(8.38)); Kk = 0.1*(1 + Cb0/(18.5));

% non-dimensional potentials
Phi_bc= vb-vc;
Phi_ca= vc-va;
Phi_lc= vl-vc*ones(1,npts); % vector
Phi_la= vl(npts)-va; 

%%% ion flux across each channel/transporter on basal and apical membrane
Jpump = Ppump*(cc0/(cc0 + K0St))^3*(cb1/(cb1 + K1St))^2;
Jn2b = Pn2b*(log(cb0*cb3^2/cc0/cc3^2)-Phi_bc); % + out of the cell
Jnkcc = Pnkcc*log(cb0*cb1*cb2^2/cc0/cc1/cc2^2);
Jae = Pae*log(cb2*cc3/cc2/cb3); 
Jkirb = C*Pkirb*(vb-vc)*(cb1 - cc1*exp(-Phi_bc))/(1 - exp(-Phi_bc));
Jn3b = Pn3b*(log(cc0*cc3^3/ca0/ca3^3)-2*Phi_ca);% + out of the cell
Jaea=Paea*log(ca2*cc3/cc2/ca3);
Jkira=C*Pkira*(Phi_ca)*(cc1 - ca1*exp(-Phi_ca))/(1 - exp(-Phi_ca));
Jcftr=-C*Pcftr*(Phi_ca)*(cc2 - ca2*exp(Phi_ca))/(1 - exp(Phi_ca));

% overall ion fluxes basal to cell through the basal membrane
J0bc=-3*Jpump+Jnkcc+Jn2b;
J1bc=2*Jpump+Jnkcc+Jkirb;
J2bc=2*Jnkcc+Jae;
J3bc=2*Jn2b-Jae;

% overall ion fluxes cell to apical through the apical membrane
J0ca=Jn3b;
J1ca=Jkira;
J2ca=Jcftr-Jaea;
J3ca=3*Jn3b+Jaea;

end

