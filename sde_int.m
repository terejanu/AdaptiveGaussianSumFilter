function X_pf = sde_int(model, time, tint, x_ini)

X_pf = x_ini;
no_particles = size(x_ini,2);
for j = 1 : time.dt_N
    X_pf = X_pf + time.dt_sde*feval(model.fx, tint(1)+time.dt_sde*(j-1), X_pf);
    X_pf = X_pf + sqrt(time.dt_sde*model.Q)*randn(model.fn,no_particles);
end;