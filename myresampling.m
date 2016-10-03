function J = resampling(q)
%Systematic Resampling or Deterministic Resampling
N = length(q);

c = cumsum(q);

J = zeros(1,N);


i = 1;
u1 = rand/N;

for j=1:N,
    u = u1 + (j-1)/N;
    while u>c(i),
        i = i + 1;
    end
    J(j) = i;
    
end
    