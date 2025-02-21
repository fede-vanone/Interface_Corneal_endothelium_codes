function [bvperr,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10] = ...
    solve_end_cleft(npts,xx,cc0,cc1,cc2,cc3,va,vc);
% this function takes in input the cellular solution and computes the
% solution of the transport problem in the cleft

global F R T Lp Ltj Ppump Pnkcc Pkira Pcftr Pkirb Pn2b Pn3b Pae Paea zx Ptj Atj P1tj P2tj P3tj...
 ca0 ca1 ca2 ca3 Aa ...
 cb0 cb1 cb2 cb3 Ab vb ...
 L h Ux Q C...
 D0 D1 D2 D3 ...
 xmesh 

%% solution initialization and settings
y_init=[cb0,0,cb1,0,cb2,0,vb,0,0,0];
solinit=bvpinit(xmesh,y_init);
options = bvpset('RelTol',1e-12,'Stats','off');

%% solution
sol=bvp5c(@(x,y)cleft_odes_coupling(x,y,cc0,cc1,cc2,cc3,va,vc,xx),...
    @(ya,yb)cleft_bc_coupling(ya,yb,cc0,cc1,cc2,cc3,va,vc,xx),solinit,options);
s=sol.stats;
bvperr=s(1).maxerr;

%% output

% concentration cl3 by electroneutrality
output_on_xmesh=deval(sol,xmesh);

y1=output_on_xmesh(1,:); % Cl0
y2=output_on_xmesh(2,:); % dCl0
y3=output_on_xmesh(3,:); % Cl1
y4=output_on_xmesh(4,:); % dCl1
y5=output_on_xmesh(5,:); % Cl2
y6=output_on_xmesh(6,:); % dCl2
y7=output_on_xmesh(7,:); % V
y8=output_on_xmesh(8,:); % dV
y9=output_on_xmesh(9,:); % p
y10=output_on_xmesh(10,:); % dp

end
