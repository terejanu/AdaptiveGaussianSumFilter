function ll = getLoglikelihood(X,mu,sig,weig)

ll = 0;

for i = 1 : size(X,2)
    s = 0;
    for j = 1 : length(weig)
        s = s + weig(j)*getLikelihood(mu{j}-X(:,i), sig{j});
    end
    if s ~= 0
        ll = ll + log(s);
    end
end