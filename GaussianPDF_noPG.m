function [pg, dpgdx, ddpgdx, dpgdmu, dpgdsig] = GaussianPDF_noPG(x,mu,sig)

r = x - mu;
pg = getLikelihood(r, sig);
isig = pinv(sig);

dpgdmu = isig*r; % <-- pg has been removed
dpgdsig = 1/2 * (isig*r*r'*isig - isig); % <-- pg has been removed
dpgdx = -dpgdmu;
ddpgdx = 2*dpgdsig;
