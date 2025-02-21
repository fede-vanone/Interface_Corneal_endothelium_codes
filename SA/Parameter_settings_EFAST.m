%% PARAMETER INITIALIZATION
% set up max and mix matrices
%%folding factors

[F,R,T,h,L,Ux,Q,C,zx,Lp,Ltj,...
    Atj,w,D0,D1,D2,D3,cb0,cb1,cb2,cb3,vb,ca0,ca1,ca2,ca3,...
    toll,nit,npts,epsilon,mu,lambdaD,sigma0,Cm]=pars5();
ref=1e-6;
% baseline parameters from the 0D optimization ( new % old)

Ppump= 4.197716594624170e-06; % 4.239302792131800e-06; 
% Pnkcc= 3.519978970519485e-06; % 3.688801964419552e-06; 
% Pkirb= 2.074799294625940e-07; % 1.751734491840590e-07; 
% Pkira= 1.544574301723440e-07; % 2.030772731366932e-07; 
% Pcftr= 1.743852885136938e-07; % 1.857571870914809e-07; 
% Pn2b= 4.661932732935699e-06; % 4.015177290133518e-06;
% Pn3b= 2.972204545461175e-06; % 2.704784522349907e-06; 
Paea= 5.318806049015642e-07; % 5.265350395316532e-07; 
Pae= 4.547240291363707e-07; % 4.937817724439357e-07;
% Ptj= 0.002438293282052; % 0.002096591900218; 
% P1tj= 0.003077082580040; % 0.002567732326705;
% P2tj= 0.001131927484887; % 0.001061168908877; 
% P3tj= 0.001002449976630; % 0.001275318160667; 

% Typical concentrations


pmean = [3.519978970519485e-06, % 3.688801964419552e-06,% Pnkcc
    2.074799294625940e-07, % 1.751734491840590e-07, % Pkirb
    1.544574301723440e-07, % 2.030772731366932e-07, % Pkira
    1.743852885136938e-07, % 1.857571870914809e-07, % Pcftr
    4.661932732935699e-06, % 4.015177290133518e-06, % Pn2b
    2.972204545461175e-06, % 2.704784522349907e-06, % Pn3b
    0.002438293282052, % 0.002096591900218, % Ptj
    0.003077082580040, % 0.002567732326705, % P1tj
    0.001131927484887, % 0.001061168908877, % P2tj
    0.001002449976630, % 0.001275318160667, % P3tj
    1]; 

pmax = pmean*1.5;
pmin = pmean*0.5; 

% Parameter Labels
efast_var={'Pnkcc','Pkirb','Pkira','Pcftr','Pn2b','Pn3b','Ptj','P1tj','P2tj','P3tj','dummy'};

% PARAMETER BASELINE VALUES (in terms of mean values)
%%% Permeabilities
% pmean=(pmax+pmin)/2;


y0=1:17;


% cc0,cc1,cc2,cc3,vc,va,xx,J0bc,J1bc,J2bc,J3bc,J0ca,J1ca,J2ca,J3ca,Jtj,J1tj,J2tj,J3tj
% Variables Labels
y_var_label={'cc0','cc1','cc2','cc3','vc','va','xx','J0tot','J1tot','J2tot',...
    'J3tot','Jtj','J1tj','J2tj','J3tj','overall water flux','water flux hexagon'};

