function [x,P] = merge(w,mu,sig)
%
% Get the mean and the covariance of the Gaussian Sum by
% merging all the Gaussian components
%
% x   - estimate the mean
% P   - estimate the covariance
%
% Gabriel Terejanu (terejanu@buffalo.edu)

%% ------------------------------------------------------------------------
% init
%--------------------------------------------------------------------------
ngs = length(w);    % number of Gaussian components
n = size(mu{1},1);  % dimensionality of the state space

x = zeros(n,1);
P = zeros(n,n);

%% ------------------------------------------------------------------------
% compute the mean and the covariance
%--------------------------------------------------------------------------
for j = 1 : ngs
    x = x + w(j)*mu{j};
end;

for j = 1 : ngs
    P = P + w(j)*(sig{j} + (mu{j}-x)*(mu{j}-x)');
end;