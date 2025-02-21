function dydx=cleft_odes_coupling(x,y,cc0,cc1,cc2,cc3,va,vc,xx)
% this function writes the cleft odes for bvp5c

global F R T Lp Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx  Atj ...
 ca0 ca1 ca2 ca3  ...
 cb0 cb1 cb2 cb3  vb ...
 L h Ux Q C...
 D0 D1 D2 D3 ...
 epsilon mu lambdaD sigma0 Cm ...
 chi bar_sigma0 H fac half_fac npts D

%% initialization

chi=Cm*lambdaD/epsilon;

bar_sigma0=sigma0*F*lambdaD/epsilon/R/T;
fac=2*Ux*L;
half_fac=fac/2;
H= h/2/h; % dimensionless y coordinate of the lateral cell membrane

% variables
y1=y(1); % c0
y2=y(2); % dc0/dx
y3=y(3); % c1
y4=y(4); % dc1/dx
y5=y(5); % c2
y6=y(6); % dc2/dx
y7=y(7); % V
y8=y(8); % dV/dx
y9=y(9); % p
y10=y(10); % dp/dx

bar_sigma=chi*(y7-vc)-bar_sigma0;

% dissociation constants for Na/K pump
K0=0.2*(1+cc1*C/8.33); % dimensional
K0St=K0/C; % dimensionless

K1=0.1*(1+y1*C/18.5); % dimensional
K1St=K1/C; % dimensionless

%% fluxes

% ion fluxes across each channel/transporter across the lateral membrane
Jpump=Ppump*cc0^3*y3^2/(K0St+cc0)^3/(K1St + y3)^2;
Jnkcc=Pnkcc*log(y1*y3*y5^2/cc0/cc1/cc2^2);
Jn2b= Pn2b*(log(y1*(y1+y3-y5)^2/cc0/cc3^2)-(y7-vc));
Jkirb=Pkirb*C*(y7-vc)*(y3-cc1*exp(vc-y7))/(1-exp(vc-y7));
Jae=Pae*log(y5*cc3/cc2/(y1+y3-y5));

% overall ion fluxes through from lateral to cell through lateral membrane
j0lc=-3*Jpump+Jnkcc+Jn2b;
j1lc=2*Jpump+Jnkcc+Jkirb;
j2lc=2*Jnkcc+Jae;
j3lc=2*Jn2b-Jae;

v_osm_lc=Lp*(R*T*C*L*(xx+cc0+cc1+cc2+cc3-2*(y1+y3))/Ux/h+mu*L^2*y9/h^3);

% slip velocity
u_slip=-bar_sigma*y8/sqrt(2*(y1+y3))-bar_sigma^2*(y2+y4)/16/(y1+y3)^2;

% second derivative of the potential from electroneutrality
d2V=(-(y2+y4)*y8-half_fac*v_osm_lc*(y1/D0+y3/D1-y5/D2-(y1+y3-y5)/D3)+...
    half_fac*H*(-H^2*y10/3+u_slip)*(y2/D0+y4/D1-y6/D2-(y2+y4-y6)/D3)+L^2*(j0lc/D0+j1lc/D1-j2lc/D2-j3lc/D3)/h/C)/(y1+y3);

% second derivative of cleft concentrations
d2c0=-y2*y8-y1*d2V+fac*(-y1*v_osm_lc-H^3*y2*y10/3+H*y2*u_slip)/D0+2*L^2*j0lc/h/C/D0;
d2c1=-y4*y8-y3*d2V+fac*(-y3*v_osm_lc-H^3*y4*y10/3+H*y4*u_slip)/D1+2*L^2*j1lc/h/C/D1;
d2c2=+y6*y8+y5*d2V+fac*(-y5*v_osm_lc-H^3*y6*y10/3+H*y6*u_slip)/D2+2*L^2*j2lc/h/C/D2;

% \hat{C} and its derivatives
CofX=2*(y1+y3);
dCofX=2*(y2+y4);
d2CofX=2*(d2c0+d2c1);

% derivative of the slip velocity
duslip=-chi.*y8.^2./sqrt(CofX)+bar_sigma.*dCofX.*y8./2./CofX./sqrt(CofX)-bar_sigma.*d2V./sqrt(CofX)...
    -bar_sigma.*chi.*y8.*dCofX./4./CofX.^2+bar_sigma.^2.*dCofX.^2./4./CofX.^3-bar_sigma.^2.*d2CofX./8./CofX.^2;

% second derivative of the pressure
d2P=v_osm_lc*3/H^3+3*duslip/H^2;

%% ODEs

dydx=[y2;d2c0;y4;d2c1;y6;d2c2;y8;d2V;y10;d2P];

end