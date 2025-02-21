function error = objectiveFunction(input_data, observed_data, Jnorm)
% this is the objective function used for the constrained ouptimization in
% 'main_in_out.m'. It takes as input the starting vector input_data, the
% reference output vector observed_data and the chosen norm (among 'lsq',
% 'mape', 'lmape', 'lmape_g'.

% run the model 
[cc0,cc1,cc2,cc3,va,vc,J0,J2,J3] = in_out(input_data);
model_output = [cc0,cc1,cc2,cc3,va,vc,J0,J2,J3];
switch Jnorm
    case 'lsq'
    error = sum((model_output(1:6) - observed_data(1:6)).^2./observed_data(1:6).^2)+0.2*sum((model_output(7:9) - observed_data(7:9)).^2./observed_data(7:9).^2) % lsq
    case 'mape'
    error = 100 * sum(  abs( (model_output - observed_data)./observed_data ))/numel(observed_data);  % mape
    case 'lmape'
    error = 100 * sum(  log(abs( (model_output - observed_data)./observed_data )))/numel(observed_data) ; % lmape
    case 'lmape_g'
    error = sum( abs( log(model_output)/log(observed_data) -1))  /numel(observed_data) ; % lmapeg
end

end