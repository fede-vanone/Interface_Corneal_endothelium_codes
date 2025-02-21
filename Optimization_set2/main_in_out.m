% This script runs the constrained optimization using fmincon. The
% objective function is written in 'objectiveFunction.m'. There are
% different options to select for the norm to use to compute the
% objective function. You should specify the reference for the input vector
% (input_ref). The optimization will use its random perturbation as a
% starting point. The values for which you want to optimize (cellular
% concentrations, cellular potential, transendothelial potential, sodium,
% chloride and bicarbonate transendothelial fluxes) are specified in
% output_ref. You can choose the algorithm used by fminsearch to find the
% optimum.

clear all
close all
%clc

rng('shuffle');

%% parameters (fixed quantities)
global  F R T Lp Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx Ptj Atj P1tj P2tj P3tj...
 ca0 ca1 ca2 ca3 Aa ...
 cb0 cb1 cb2 cb3 vb Ab  ...
 L h Ux Q C...
 xmesh w...
 D0 D1 D2 D3 ...
 Ltj w_cheb_01 ...
 ref epsilon mu lambdaD sigma0 Cm chi bar_sigma0 H fac half_fac npts D nit toll

[F,R,T,h,L,Ux,Q,C,zx,Lp,Ltj,Atj,w,D0,D1,D2,D3,cb0,cb1,cb2,cb3,vb,ca0,ca1,ca2,ca3,toll,nit,npts,epsilon,mu,lambdaD,sigma0,Cm]=pars5();

chi=Cm*lambdaD/epsilon;

bar_sigma0=sigma0*F*lambdaD/epsilon/R/T;
fac=2*Ux*L;
ref=1e-6;

% we fix the pump from set 1
Ppump_nd=4.197716594624170e-06/ref;
Ppump=Ppump_nd*ref;

%% Define ground truth
% Important: one should remove the 1e6 scaling when defining the inputs to be calibrated. 
% The /ref is done later in the in_out routine

% set 1 as input reference 
input_ref=[3.519978970519485e-06; % nkcc
    C*2.074799294625940e-07; % kirb
C* 1.544574301723440e-07; % kira
C* 1.743852885136938e-07; % cftr
4.661932732935699e-06; % n2b
2.972204545461175e-06; % n3b
5.318806049015642e-07; % aea
4.547240291363707e-07; % ae
C* 0.002438293282052; % Ptj
C* 0.003077082580040; % P1tj
C* 0.001131927484887; % P2tj
C* 0.001002449976630 % P3tj
];
% if you want to test the optimizer, use the computed output from input_ref 
% as "experimental values" for the output, and you should get back input_ref 
% from the optimization:
% [cc0_exp,cc1_exp,cc2_exp,cc3_exp,va_exp,vc_exp,J0_exp,J2_exp,J3_exp]=in_out(input_ref); 

% if you want to optimize against data, define the data here:
cc0_exp=15/C; 
cc1_exp=132/C;
cc2_exp=38/C;
cc3_exp=25/C; 
va_exp=(-.5e-3)*F/R/T;
vc_exp=(-60e-3)*F/R/T;
J0_exp=(9e-6)/ref;
% J1_exp=(2.5e-5)/ref;
J2_exp=(6e-6)/ref;
J3_exp=(5e-6)/ref;

output_ref = [cc0_exp,cc1_exp,cc2_exp,cc3_exp,va_exp,vc_exp,J0_exp,J2_exp,J3_exp];

%% Choose the initial condition by perturbing input_ref
% ->0.95 test case
% input0=input_ref*0.95;
% prange=0.07;

%->random perturbation (multiplicative number uniform in 1 +[-prange,prange])
% prange=0.5; % try different values of prange
% perturbation=(rand(numel(input_ref),1)-0.5)/0.5*prange;
% input0=input_ref.*(1 + perturbation  );

x0=input_ref/5;
x1=input_ref*5;
% starting point
input0=x0+rand(length(input_ref),1).*(x1-x0);

% if you want to start exactly from input_ref: 
% input0=input_ref;
%% Calibration
observed_data = output_ref;
%parpool
% choose the norm to use in objectiveFunction
Jnorm1='lsq';
Jnorm2='mape';
Jnorm3='lmape';
Jnorm4='lmape_g';
Jnorm=Jnorm1;
% value of the objective function with the permeabilities in input_ref and 
% the experimental data in observed_data, computed with Jnorm
J = objectiveFunction(input_ref, observed_data,Jnorm);

% algorithms for the optimization
alg1='interior-point'; 
alg3='sqp'; 
alg4='sqp-legacy'; 
alg5='active-set';
% choose the algorithm
alg=alg1;

options = optimoptions('fmincon','Display', 'iter', 'Algorithm', alg,'TolFun', 1e-3, 'TolX', 1e-6 );
% run the optimization with fmincon
[input_opt, fval] = fmincon(@(input) objectiveFunction(input, observed_data,Jnorm), ...
    input0, [], [], [], [], input_ref/5, input_ref*5, [], options);

%% Checks
% value of the objective function at the last iteration
J = objectiveFunction(input_opt, observed_data,Jnorm)
% if you used the computed output as reference, check this
% error_input=sum( abs( (input_opt - input_ref)./input_ref ) )*100;
% fprintf('J=%e   Err_inp=%e\n',J,error_input)

% look at the optimum vs the reference input
[input_opt(:),input_ref(:)]

% save the workspace, the value of the objective function at the last
% iteration and the vector for the optimum
save('optim_workspace.mat')
save('objective.mat','J')
save('optimum_perm.mat','input_opt')