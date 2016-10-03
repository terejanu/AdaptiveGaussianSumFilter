function [omu,osig,oll] = getGSdata(iweig,imu,isig,pf_samples,pf_w)
%
% Get the estimates for one particular Gaussian Sum
%
% omu   - estimate the mean
% osig  - estimate the covariance
% oll   - log-likelihood of the particles
%
% Gabriel Terejanu (terejanu@buffalo.edu)

%% ------------------------------------------------------------------------
% omu - compute the mean of the Gaussian Sum
% osig - compute the covariance of the Gaussian Sum
%--------------------------------------------------------------------------
[omu,osig] = merge(iweig,imu,isig);

%% ------------------------------------------------------------------------
% oll - compute log-likelihood of the particles
%--------------------------------------------------------------------------
no_particles = length(pf_w);

% resampling
I = myresampling(pf_w);
I = round(I);
tmp_X_pf = zeros(size(pf_samples));
for j = 1 : no_particles
    tmp_X_pf(:,j) = pf_samples(:,I(j));
end;    
    
oll = getLoglikelihood(tmp_X_pf, imu, isig, iweig);