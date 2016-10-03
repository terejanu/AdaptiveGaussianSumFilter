function [omu,osig] = getPFdata(pf_samples,pf_w)
%
% Get the estimates for the Particle Filter
%
% omu   - estimate the mean
% osig  - estimate the covariance
%
% Gabriel Terejanu (terejanu@buffalo.edu)

%% ------------------------------------------------------------------------
% init & resampling
%--------------------------------------------------------------------------
no_particles = length(pf_w);

% resampling
I = myresampling(pf_w);
I = round(I);
tmp_X_pf = zeros(size(pf_samples));
for j = 1 : no_particles
    tmp_X_pf(:,j) = pf_samples(:,I(j));
end;    

%% ------------------------------------------------------------------------
% omu - compute the mean of the Particle Filter
%--------------------------------------------------------------------------
omu = mean(tmp_X_pf,2);

%% ------------------------------------------------------------------------
% osig - compute the covariance of the Particle Filter
%--------------------------------------------------------------------------
osig = cov(tmp_X_pf');