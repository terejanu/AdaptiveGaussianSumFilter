function Jx = fx_lorenz_jac(t, x)

r1 = 10; r2 = 28; r3 = 8/3;
x1 = x(1); x2 = x(2); x3 = x(3);

Jx = [-r1, r1, 0; r2-x3, -1, -x1; x2, x1, -r3];