% this code runs the coupled model cell-cleft

global C ca0 ca1 ca2 ca3 cb0 cb1 cb2 cb3 vb npts H chi bar_sigma0 
%% Parameters

[F,R,T,h,L,Ux,Q,C,zx,Lp,Ltj,...
    Atj,w,D0,D1,D2,D3,cb0,cb1,cb2,cb3,vb,ca0,ca1,ca2,ca3,...
    toll,nit,npts,epsilon,mu,lambdaD,sigma0,Cm]=pars5();
% Permeabilities are taken from the X matrix

ref=1e-6;

% Peclet number 
Pe=L*Ux/D0;
fprintf('Peclet number=%g\n',Pe)

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

%% new version: start with one cell solution with the variables initialized asa in the basal
[cc0_start,cc1_start,cc2_start,cc3_start,vc_start,va_start,xx_start,J0bc,J1bc,J2bc,J3bc,J0ca,J1ca,J2ca,J3ca]=...
    cell_full(npts,cl0,cl1,cl2,cl3,vl,p);

cc0=cc0_start; cc1=cc1_start; cc2=cc2_start; cc3=cc3_start; xx=xx_start; ...
    va=va_start; vc=vc_start;

%% iteration 
tic
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

% compute te errors
% relative cell error 
err_cell_vec=[abs(cc0-cc0_n)/cc0,abs(cc1-cc1_n)/cc1, abs(cc2-cc2_n)/cc2, abs(cc3-cc3_n)/cc3,...
    abs(va-va_n)/abs(va), abs(vc-vc_n)/abs(vc), abs(xx-xx_n)/xx];
err_cell_norm(n)=max(err_cell_vec); % maximum norm
% err_cell_norm(n)=vecnorm(err_cell_vec); % 2-norm

% relative cleft error: each row contains the error on the concentration of
% 1 solute/the potential/ pressure
err_cleft_mat=[abs(cl0-cl0_n)/vecnorm(cl0);abs(cl1-cl1_n)/vecnorm(cl1);abs(cl2-cl2_n)/vecnorm(cl2);...
    abs(cl3-cl3_n)/vecnorm(cl3);abs(vl-vl_n)/vecnorm(vl);abs(p_n-p)/vecnorm(p)]; % 5*npts
 
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
cc0=cc0_n;cc1=cc1_n; cc2=cc2_n; cc3=cc3_n; vc=vc_n; xx=xx_n; va=va_n; %qx=qx_n;

% use directly these fluxes, because they don't depend on the coupling
J0bc=J0bc_n; J1bc=J1bc_n; J2bc=J2bc_n; J3bc=J3bc_n; 
J0ca=J0ca_n; J1ca=J1ca_n; J2ca=J2ca_n; J3ca=J3ca_n; 

% non-dimensional potentials
Phi_ca= vc-va;
Phi_lc= vl-vc*ones(1,npts); % vector
Phi_la= vl(npts)-va; 

% tight junction fluxes
Jtj=Atj*Ptj*(Phi_la)*C*(cl0(npts) - ca0*exp(-Phi_la))/(1-exp(-Phi_la));
J1tj=Atj*P1tj*(Phi_la)*C*(cl1(npts) - ca1*exp(-Phi_la))/(1-exp(-Phi_la));
J2tj=-Atj*P2tj*(Phi_la)*C*(cl2(npts) - ca2*exp(Phi_la))/(1-exp(Phi_la));
J3tj=-Atj*P3tj*(Phi_la)*C*(cl3(npts) - ca3*exp(Phi_la))/(1-exp(Phi_la));


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

bar_sigma=chi*(y7-vc)-bar_sigma0;

% water velocity from lateral to cell, osmotic + hydraulic pressure
v_osm_lc=Lp*(R*T*C*L*(xx+cc0+cc1+cc2+cc3-2*(y1+y3))/Ux/h+mu*L^2*y9/h^3);

% \hat{C} and its derivative
CofX=2*(y1+y3);
dCofX=2*(y2+y4);
% slip velocity
u_slip=-bar_sigma.*y8./sqrt(CofX)-bar_sigma.^2.*dCofX./8./CofX.^2;

%% fluxes and conservations

% dimensional concentrations
Cb0=C*cb0; Cb1=C*cb1; Cb2=C*cb2; Cb3=C*cb3; 
Cc0=C*cc0; Cc1=C*cc1; Cc2=C*cc2; Cc3=C*cc3; 
Ca0=ca0*C; Ca1=ca1*C; Ca2=ca2*C; Ca3=ca3*C; Xx=xx*C;

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
Qlat=qlat*Ux*h; % dimensional m^2/s
overall_water_flux=(qtj*Ux*h+Qca)/(w+h); % dimensional water flux towards the apical compartment, in m/s

J0tot=J0ca+Jtj;
J1tot=J1ca+J1tj;
J2tot=J2ca+J2tj;
J3tot=J3ca+J3tj;

%% fluxes including area factors

Acell=(w/2)^2*3*sqrt(3)/2;
Atot=(w/2+h/sqrt(3))^2*3*sqrt(3)/2;
Acleft=Atot-Acell;

uTJ=Qtj/h; % velocity of water at the tight junction
uca=Qca/w; % velocity of water at the apical cell membrane
% mean velocity across the endothelium accounting for the area factors
Qcell=uca*Acell;
Qcleft=uTJ*Acleft;
u_end_hex=1/Atot*(Qcell+Qcleft);