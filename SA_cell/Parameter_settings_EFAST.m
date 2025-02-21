%% PARAMETER INITIALIZATION
% set up max and mix matrices
%%folding factors

[F,R,T,h,L,Ux,Q,C,zx,Lp,Ltj,...
    Atj,w,D0,D1,D2,D3,cb0,cb1,cb2,cb3,vb,ca0,ca1,ca2,ca3,...
    toll,nit,npts,epsilon,mu,lambdaD,sigma0,Cm]=pars5();
ref=1e-6;

% baseline parameters from the set 1
% fix the pump
Ppump= 4.197716594624170e-06;
% define the mean values for the parameter distributions
pmean = [3.519978970519485e-06, 
    2.074799294625940e-07, 
    1.544574301723440e-07, 
    1.743852885136938e-07, 
    4.661932732935699e-06, 
    2.972204545461175e-06, 
    5.318806049015642e-07,
    4.547240291363707e-07, 
    0.002438293282052, 
    0.003077082580040, 
    0.001131927484887, 
    0.001002449976630, 
    1]; 

% ranges for the uniform distributions
pmax = pmean*1.5;
pmin = pmean*0.5; 

% Parameter Labels
efast_var={'Pnkcc','Pkirb','Pkira','Pcftr','Pn2b','Pn3b','Paea','Pae','Ptj','P1tj','P2tj','P3tj','dummy'};

% length of the vector of variables
y0=1:19;

% Variables Labels
y_var_label={'cc0','cc1','cc2','cc3','vc','va','xx','J0bc','J1bc','J2bc',...
    'J3bc','J0ca','J1ca','J2ca','J3ca','Jtj','J1tj','J2tj','J3tj'};

