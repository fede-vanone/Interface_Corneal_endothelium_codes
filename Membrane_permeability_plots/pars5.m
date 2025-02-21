%% PARAMETER INITIALIZATION
function [F,R,T,h,L,Ux,Q,C,zx,...
    Atj,w,D0,D1,D2,D3,cb0,cb1,cb2,cb3,vb,ca0,ca1,ca2,ca3,...
    toll,nit,npts,epsilon,mu,lambdaD]=pars5()
% this function assigns the values to the parameters

%% thermodynamical constants
F = 96485.332; 
R = 8.314;
T = 34 + 273;
epsilon=73*8.854e-12; % dielectric permittivity
mu= 7.34e-4; % viscosity % viscosity for water at 34°C (Pa s)

%% cleft dimensions
h=30e-9; % width
L=5e-6; % length

%% width of the cell
w=20e-6;

%% scales
Ux=epsilon*R^2*T^2/L/mu/F^2;
Q=Ux*h; % flux m^2/s
C=100; % concentration mol/m^3

lambdaD=sqrt(epsilon*R*T/F^2/C); % Debye length

%% Valence of the fixed charges in the cell
zx=-1.05;

%% hydraulic conductivity of the membrane specified in 'plots_Lp.m' or 'plots_Ltj.m'

%% tight junction permeability to water specified in 'plots_Lp.m' or 'plots_Ltj.m'

%% area factor for the tj
Atj=h/w;

%% diffusion coefficients
D0=3e-9; D1=2.22e-9; D2=1.69e-9; D3=1.36e-9; 

%% basal dimensionless concentrations and potential
cb0= 140.1/C; % Na
cb1 = 4.9/C; % K
cb3= 37/C; % HCO3
%Cb2 = Ca0+Ca1-Ca3; % Cl
cb2 = 108/C;

% basal potential
vb=0*F/R/T;


%% apical concentrations
ca0= 140.1/C; % Na
ca1 = 4.9/C; % K
ca3= 37/C; % HCO3
%Ca2 = Ca0+Ca1-Ca3; % Cl
ca2 = 108/C;

%% integration
% tolerance (reduced to speed up the computations)
toll=1e-6;

% max number of iterations 
nit=150;

% number of mesh points
npts=100; 
end