function [cc0,cc1,cc2,cc3,va,vc,J0,J2,J3]=in_out(input)
% This function takes a vector of 12 permeabilities as input and gives as
% output the cellular concentrations, cellular and transendothelial
% potential, sodium, chloride and bicarbonate transendothelial fluxes computed by the model.

global  F R T Lp Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx Ptj Atj P1tj P2tj P3tj...
 ca0 ca1 ca2 ca3 Aa ...
 cb0 cb1 cb2 cb3 vb Ab  ...
 L h Ux Q C...
 xmesh w...
 D0 D1 D2 D3 ...
 Ltj w_cheb_01 ...
 ref epsilon mu lambdaD sigma0 Cm chi bar_sigma0 H fac half_fac npts D nit toll

%%
% make the permeablities dimensionless
Pnkcc_nd=input(1)/ref;
Pkirb_nd=input(2)/ref;
Pkira_nd=input(3)/ref;
Pcftr_nd=input(4)/ref;
Pn2b_nd=input(5)/ref;
Pn3b_nd=input(6)/ref;
Paea_nd=input(7)/ref;
Pae_nd=input(8)/ref;
Ptj_nd=input(9)/ref;
P1tj_nd=input(10)/ref;
P2tj_nd=input(11)/ref;
P3tj_nd=input(12)/ref;


% dimensional permeabilities 
Pnkcc=Pnkcc_nd*ref;
Pkirb=Pkirb_nd*ref/C;
Pkira=Pkira_nd*ref/C;
Pcftr=Pcftr_nd*ref/C;
Pn2b=Pn2b_nd*ref;
Pn3b=Pn3b_nd*ref;
Paea=Paea_nd*ref;
Pae=Pae_nd*ref;
Ptj=Ptj_nd*ref/C;
P1tj=P1tj_nd*ref/C;
P2tj=P2tj_nd*ref/C;
P3tj=P3tj_nd*ref/C;

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
% p2 per non aver bisgno del pde toolbox
[w_cheb,x_cheb]=clenshaw_curtis_p2(npts,1); % x_cheb in [1,-1]
xmesh=flip(x_cheb+1)'/2; % flip [-1,1], shift to [0,2] and rescale to [0,1]
mapping_01=1/2;
w_cheb_01=flip(w_cheb')*mapping_01; % weights for integration formula in 01
[x_chebdif, DM] = chebdif(npts, 1); % x_chebdif same as x_cheb (up to 10^-16)
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

% start with 1 cell iteration
[cc0_start,cc1_start,cc2_start,cc3_start,vc_start,va_start,xx_start,J0bc,J1bc,J2bc,J3bc,J0ca,J1ca,J2ca,J3ca]=...
    cell_full(npts,cl0,cl1,cl2,cl3,vl,p);

cc0=cc0_start; cc1=cc1_start; cc2=cc2_start; cc3=cc3_start; xx=xx_start; ...
    va=va_start; vc=vc_start;

%% iteration 

for n=1:nit
  % fprintf('Iteration %g \n',n) 

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

% solve for the cell
[cc0_n, cc1_n,cc2_n, cc3_n, vc_n, va_n, xx_n,J0bc_n,J1bc_n,J2bc_n,J3bc_n,J0ca_n,J1ca_n,J2ca_n,J3ca_n]...
    = cell_full(npts,cl0_n,cl1_n,cl2_n,cl3_n,vl_n,p_n);

% relative cell error 
err_cell_vec=[abs(cc0-cc0_n)/cc0,abs(cc1-cc1_n)/cc1, abs(cc2-cc2_n)/cc2, abs(cc3-cc3_n)/cc3,...
    abs(va-va_n)/abs(va), abs(vc-vc_n)/abs(vc), abs(xx-xx_n)/xx];
err_cell_norm(n)=max(err_cell_vec); % maximum norm
% err_cell_norm(n)=vecnorm(err_cell_vec); % 2-norm

% relative cleft error: each row contains the error on the concentration of
% 1 solute/the potential/ the pressure
err_cleft_mat=[abs(cl0-cl0_n)/vecnorm(cl0);abs(cl1-cl1_n)/vecnorm(cl1);abs(cl2-cl2_n)/vecnorm(cl2);...
    abs(cl3-cl3_n)/vecnorm(cl3);abs(vl-vl_n)/vecnorm(vl);abs(p_n-p)/vecnorm(p)]; % 5*npts
% err_cleft_mat=[(Cl0-Cl0_n)/C;(Cl1-Cl1_n)/C;(Cl2-Cl2_n)/C;...
%    (Cl3-Cl3_n)/C;(Vl-Vl_n)*F/R/T];
 
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
y1=y1_n;
y2=y2_n;
y3=y3_n;
y4=y4_n;
y5=y5_n;
y6=y6_n;
y7=y7_n;
y8=y8_n;
y9=y9_n;
y10=y10_n;

% cell
cc0=cc0_n;cc1=cc1_n; cc2=cc2_n; cc3=cc3_n; vc=vc_n; xx=xx_n; va=va_n; 

% use directly these fluxes, because they don't depend on the coupling
J0bc=J0bc_n; J1bc=J1bc_n; J2bc=J2bc_n; J3bc=J3bc_n; 
J0ca=J0ca_n; J1ca=J1ca_n; J2ca=J2ca_n; J3ca=J3ca_n; 

%%% check equations for the cell 
% (exactly those solved, because after the cell the varibles were not
% updated again)

% non-dimensional potentials
% Phi_bc= vb-vc;
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

% 
Jpump_l = Ppump*(cc0/(cc0 + K0st))^3*(cl1./(cl1 + K1_l_st)).^2;
Jn2b_l = Pn2b*(log(cl0.*cl3.^2/cc0/cc3^2)-Phi_lc); % + out of the cell
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
% uy=Lp*R*T*C*(xx+cc0+cc1+cc2+cc3-cl0-cl1-cl2-cl3);
uy_nd=xx+cc0+cc1+cc2+cc3-cl0-cl1-cl2-cl3;
% eq_wb = (xx+cc0+cc1+cc2+cc3-cb0-cb1-cb2-cb3)+2*L*sum(w_cheb_01.*uy_nd)/w-...
%     (ca0+ca1+ca2+ca3-cc0-cc1-cc2-cc3-xx);
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

%%% stopping criteria
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
% compute the transendothelial ion fluxes
J0=(Jtj+J0ca)/ref; J1=(J1tj+J1ca)/ref; J2=(J2tj+J2ca)/ref; J3=(J3tj+J3ca)/ref;


