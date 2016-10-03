function y = fx_lorenz(t, x)

r1 = 10; r2 = 28; r3 = 8/3;
x1 = x(1,:); x2 = x(2,:); x3 = x(3,:);

y = zeros(size(x));

y(1,:) = r1*(-x1 + x2);
y(2,:) = r2*x1 - x2 - x1.*x3;
y(3,:) = -r3*x3 + x1.*x2;