function sigma = IUT_getSigmaPoints(x, S, n)

alpha = normpdf(2)/normpdf(1);
beta = normpdf(3)/normpdf(1);

% # of sigma points to be generated = 6*n+1

sigma(1).x = x;
sigma(1).W = 1 - n*(1+alpha+beta)/(1+4*alpha+9*beta);

for i = 1 : n
    sigma(i+1).x = x + S(:,i);
    sigma(i+1).W = 1/(2+8*alpha+18*beta);
    sigma(i+n+1).x = x - S(:,i);
    sigma(i+n+1).W = 1/(2+8*alpha+18*beta);

    sigma(i+2*n+1).x = x + 2*S(:,i);
    sigma(i+2*n+1).W = alpha/(2+8*alpha+18*beta);
    sigma(i+3*n+1).x = x - 2*S(:,i);
    sigma(i+3*n+1).W = alpha/(2+8*alpha+18*beta);    

    sigma(i+4*n+1).x = x + 3*S(:,i);
    sigma(i+4*n+1).W = beta/(2+8*alpha+18*beta);
    sigma(i+5*n+1).x = x - 3*S(:,i);
    sigma(i+5*n+1).W = beta/(2+8*alpha+18*beta);      
end
