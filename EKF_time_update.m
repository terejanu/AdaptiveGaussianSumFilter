function [omu, osig] = EKF_time_update(myfilter, model, imu, isig, T, opt)
%
% EKF Time Update
%
% Gabriel Terejanu (terejanu@buffalo.edu)

n = model.fn;

myX = [imu; reshape(isig, numel(isig), 1)];
[t,y] = ode45(@EKF_propagation, T, myX, opt, model);

myY = y(end,:);
omu = myY(1:n)';
osig = reshape(myY(n+1:end),n,n);