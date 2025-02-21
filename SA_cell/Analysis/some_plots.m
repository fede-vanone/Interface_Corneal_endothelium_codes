% for the output corresponding to x=l, see its correlation with each of the 
% 13 parameters 

clear all
close all

load('X.mat');
load('Yfull.mat','Y');

%%
% check sampling of parameter k vs parameter m
k=5;
m=9;

% fourth index in X is number od resamplings
figure
plot(X(:,k,k,1),X(:,m,k,1),'o')
hold on 
plot(X(:,k,k,2),X(:,m,k,2),'o')
xlabel(efast_var(k))
ylabel(efast_var(m))

%%
for k=1:13 % number of parameters (permeabilities)
l=1; % which output we are looking at among the following:
% 1'cc0',2'cc1',3'cc2',4'cc3',5'vc',6'va',7'xx',8'J0bc',9'J1bc',10'J2bc',
%    11'J3bc',12'J0ca',13'J1ca',14'J2ca',15'J3ca',16'Jtj',17'J1tj',18'J2tj',19'J3tj';

% for 3 sampling curves
ind1=isnan(Y(:,1,k,1));
ind2=isnan(Y(:,1,k,2));
ind3=isnan(Y(:,1,k,3));

% fit with a 2nd order poly X,Y of the 3 resamplings
p = polyfit([X(~ind1,k,k,1);X(~ind2,k,k,2);X(~ind3,k,k,3)],[Y(~ind1,l,k,1);Y(~ind2,l,k,2);Y(~ind3,l,k,3)],2);

figure
plot(X(:,k,k,1),Y(:,l,k,1),'o')
hold on 
plot(X(:,k,k,2),Y(:,l,k,2),'o') 
plot(X(:,k,k,3),Y(:,l,k,3),'o') 
% plot the polynomial
plot(sort(X(:,k,k,1)),polyval(p,sort(X(:,k,k,1))),'k-') 
set(gca, 'FontSize', 16)
xlabel( efast_var(k)+ " (m s$^{-1}$)", 'FontSize',16,'Interpreter','LaTex')
ylabel(y_var_label(l),'FontSize',16,'Interpreter','LaTex')
end
