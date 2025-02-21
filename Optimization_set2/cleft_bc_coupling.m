function res=cleft_bc_coupling(ya,yb,cc0,cc1,cc2,cc3,va,vc,xx)
% this function sets the boundary conditions for the cleft odes

global F R T Lp Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx Ptj Atj P1tj P2tj P3tj Ltj...
 ca0 ca1 ca2 ca3...
 cb0 cb1 cb2 cb3 vb ...
 L h Ux Q C...
 D0 D1 D2 D3 ...
 epsilon mu lambdaD sigma0 Cm...
 chi bar_sigma0 H Atj fac half_fac npts

%% initialization

% variables for the basal (x=0) and apical (x=1) side of the cleft
ya1=ya(1); % c0(x=0)
ya2=ya(2); % dc0/dx (x=0)
ya3=ya(3); % c1(x=0)
ya4=ya(4); % dc1/dx (x=0)
ya5=ya(5); % c2(x=0)
ya6=ya(6); % dc2/dx (x=0)
ya7=ya(7); % V(x=0)
ya8=ya(8); % dV/dx (x=0)
ya9=ya(9); % p (x=0)
ya10=ya(10); % dp (x=0)

yb1=yb(1); % c0(x=1)
yb2=yb(2); % dc0/dx (x=1)
yb3=yb(3); % c1(x=1)
yb4=yb(4); % dc1/dx (x=1)
yb5=yb(5); % c2(x=1)
yb6=yb(6); % dc2/dx (x=1)
yb7=yb(7); % V(x=1)
yb8=yb(8); % dV/dx (x=1)
yb9=yb(9); % p (x=1)
yb10=yb(10); % dp (x=1)

bar_sigma=chi*(yb7-vc)-bar_sigma0;

u_slip_b=-bar_sigma*yb8/sqrt(2*(yb1+yb3))-bar_sigma^2*(yb2+yb4)/16/(yb1+yb3)^2;

v_osm_tj=Ltj*(R*T*C*(ca0+ca1+ca2+ca3-2*(yb1+yb3))/Ux+mu*L*yb9/h^2);

% no Atj because bc are pointwise
J0tj=Ptj*C*(yb7-va)*(yb1-ca0*exp(va-yb7))/(1-exp(va-yb7));
J1tj=P1tj*C*(yb7-va)*(yb3-ca1*exp(va-yb7))/(1-exp(va-yb7));
J2tj=-P2tj*C*(yb7-va)*(yb5-ca2*exp(yb7-va))/(1-exp(yb7-va));
J3tj=-P3tj*C*(yb7-va)*(yb1+yb3-yb5-ca3*exp(yb7-va))/(1-exp(yb7-va));

%% boundary conditions 

res=[ya1-cb0;...
    ya3-cb1;... % check the fluxes across the tight junction
    ya5-cb2;...
    -yb2-yb1*yb8+half_fac*yb1*(-H^2*yb10/3+u_slip_b)/D0-J0tj*L/D0/C;...
    -yb4-yb3*yb8+half_fac*yb3*(-H^2*yb10/3+u_slip_b)/D1-J1tj*L/D1/C;...
    -yb6+yb5*yb8+half_fac*yb5*(-H^2*yb10/3+u_slip_b)/D2-J2tj*L/D2/C;...
    -(yb2+yb4-yb6)+(yb1+yb3-yb5)*yb8+half_fac*(yb1+yb3-yb5)*(-H^2*yb10/3+u_slip_b)-J3tj*L/D3/C;
    ya7;
    ya9;
    yb10-3*(u_slip_b-v_osm_tj)/H^2];

end