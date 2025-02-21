% First order and total effect indices for a given
% model computed with Extended Fourier Amplitude
% Sensitivity Test (EFAST).
% Andrea Saltelli, Stefano Tarantola and Karen Chan.
% 1999. % "A quantitative model-independent method for global
% sensitivity analysis of model output". % Technometrics 41:39-56
clear;
close all;
rng("twister")
%% INPUT
NR = 3; %: no. of search curves - RESAMPLING
k = 10 + 1; % # of input factors (parameters varied) + dummy parameter
NS = 1000; % # of samples per search curve >=500
wantedN=NS*k*NR; % wanted no. of sample points

% OUTPUT
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of k frequencies
% OMi : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output

MI =5; %: maximum number of fourier coefficients
% that may be retained in calculating the partial
% variances without interferences between the
% assigned frequencies

%% PARAMETERS AND ODE SETTINGS (they are included in the following file)
Parameter_settings_EFAST;

% Computation of the frequency for the group
% of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
    '65 per factor.\n']);
    return;
end


%% Pre-allocation of the output matrix Y
%% Y will save only the points of interest 
Y(NS,length(y0),length(pmin),NR)=0;  % pre-allocation

% Loop over k parameters (input factors)
for i=1:k % i=# of replications (or blocks)
    % Algorithm for selecting the set of frequencies.
    % OMci(i), i=1:k-1, contains the set of frequencies
    % to be used by the complementary group.
    OMci = SETFREQ(k,OMi/2/MI,i);   
    % Loop over the NR search curves.
    for Ll=1:NR
        % Setting the vector of frequencies OM
        % for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar
        % variable S and the coordinates
        % {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,Ll) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1
        
        % Transform distributions from standard
        % uniform to general.
        X(:,:,i,Ll) = parameterdist(X(:,:,i,Ll),pmax,pmin,0,1,NS,'unif'); %% this is what assigns 'our' values rather than 0:1 dist
        
    end 
end
save X.mat