function l = getLikelihood(r, A)

l = exp(-r'*inv(A)*r/2)/sqrt(det(2.*pi.*A));