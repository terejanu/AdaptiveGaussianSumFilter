function [MM, NN] = ComputeFPE_QP_UT_noPG(weig, mu, sig, model, time, tk)

% compute matrix MM
MM = zeros(length(weig));
MMtmp = MM;
for i = 1 : length(weig)
    for j = i+1 : length(weig)
        MMtmp(i,j) = getLikelihood(mu{i}-mu{j}, sig{i}+sig{j});
    end
    MM(i,i) = 1/sqrt(det(4.*pi.*sig{i}));
end
MM = MM + MMtmp + MMtmp';  
MM = MM/time.dt^2;

% compute matrix NN
NN = zeros(length(weig));    

for i = 1 : length(weig)
    for j = 1 : length(weig)

        SigC = inv(inv(sig{i}) + inv(sig{j}));
        muC = SigC*(inv(sig{i})*mu{i} + inv(sig{j})*mu{j});
        cC = getLikelihood(mu{i}-mu{j},sig{i}+sig{j});
        
%         sigma = SRUKF_getSigmaPoints(muC, chol(SigC)', model.fn, 10^-4);
        sigma = IUT_getSigmaPoints(muC, chol(SigC)', model.fn);

         for k = 1 : length(sigma)
            x = sigma(k).x;

            [pg, dpgdx, ddpgdx, dpgdmu, dpgdsig] = GaussianPDF_noPG(x, mu{j}, sig{j});

            myX = [mu{j}; reshape(sig{j}, numel(sig{j}), 1)];
            myDX = EKF_propagation(time.tspan(tk), myX, model);
            mudot = myDX(1:model.fn);
            sigdot = reshape(myDX(model.fn+1:end), model.fn, model.fn);

            dpgidt = dpgdmu'*mudot + trace(dpgdsig*sigdot);

            fx = feval(model.fx, time.tspan(tk), x);
            Jx = feval([model.fx '_jac'], time.tspan(tk), x);
            dopgi_dot = - dpgdx' * fx - trace(Jx) + 1/2*trace(model.Q * ddpgdx);     % pg*trace(Jx)

            Pi = 1/time.dt;                             %         Pk(j) = pg/time.dt;
            Lj = dpgidt - dopgi_dot - 1/time.dt;       %         Lk(j) = dpgidt - dopgi_dot - pg/time.dt;

            NN(i,j) = NN(i,j) + sigma(k).W*Pi*Lj;
        end;
        NN(i,j) = NN(i,j)*cC;
    end;
end;