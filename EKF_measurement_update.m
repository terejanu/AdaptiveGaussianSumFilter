function [omu, osig, z_mu, z_sig] = EKF_measurement_update(myfilter, model, imu, isig, T, opt, ymeas)
%
% EKF Measurement Update
%
% Gabriel Terejanu (terejanu@buffalo.edu)

z_mu = feval(model.hx, imu);
h = feval([model.hx '_jac'], imu);
z_sig = h*isig*h';

gain = isig * h' * inv(z_sig + model.R);
osig = isig - gain*h*isig;
omu = imu + gain*(ymeas - z_mu);
