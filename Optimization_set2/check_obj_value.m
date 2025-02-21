% This script plots the value of the objective function at the last iteration 
% saved as 'objective.mat' in each of the 20 subfolders run-i, with i from
% 1 to 20. This plot is useful to assess the convergence of the different
% optimization runs.
close all
clear all

for i=[1 2 3 4 5 6 7 8 9 10 11 12 14 15 17 18 19 20]
    load(strcat(pwd,'/run-',num2str(i),'/objective.mat'))
    plot(i,J,'ko')
    hold on
end