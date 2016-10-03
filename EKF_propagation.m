function dx = EKF_propagation(t, x, model)

n = model.fn;

mu = x(1:n);
sig = reshape(x(n+1:end),n,n);

F = feval([model.fx '_jac'], t, mu);
sig = F*sig + sig*F' + model.Q;
mu = feval(model.fx, t, mu);

dx = [mu; reshape(sig, numel(sig), 1)];