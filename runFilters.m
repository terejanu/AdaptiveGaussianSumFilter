function [xt,F1,F2,F4,PF] = runFilters(iRun,model,prior,time,measAvail,pf,inc,myfilter)
%
% Run all the 4 filters for one particular MC run. Filters:
%   1. no update when measurements are available
%   2. Bayesian classical update
%   4. Continuous update based on FPKE
%
% for each Fx/PF there will be saved:
%   Fx.x{k}                     - estimate first moment (mean) - over time
%   Fx.P{k}                     - estimate second moment (covariance) - over time
%   Fx.w{k}                     - weights for gaussian components - over time
%   Fx.ll{k}                    - log-likelihood of the particles - over time
%   Fx.pdf{k}.dim(j).mmX        - pdf for dim j (1:model.fn) - X axis
%   Fx.pdf{k}.dim(j).mmY        - pdf for dim j (1:model.fn) - Y axis
%   Fx.pdf{k}.dim(j).mmZ        - pdf for dim j (1:model.fn) - Z axis
%
% Gabriel Terejanu (terejanu@buffalo.edu)

%% ------------------------------------------------------------------------
% init
%--------------------------------------------------------------------------
path(path, genpath(inc));
opt = odeset('reltol',1e-8,'abstol',1e-8);

fprintf('- Run # %d\n',iRun);
      
%% ------------------------------------------------------------------------
% create truth & measurement
%--------------------------------------------------------------------------
fprintf('  - create truth & measurements\n');
xt = cell(1, time.nSteps);
y_meas = cell(1, time.nSteps);

% select a gaussian component
rand('state',iRun*100);
gs_sel = 1;
u = rand;
u_total = 0;
for j = 1 : prior.n
    if ((u >= u_total) && (u < u_total + prior.weig(j)))
        gs_sel = j;
        break;
    end
    u_total = u_total + prior.weig(j);
end

% draw a sample from the chosen gaussian component
randn('state',iRun*100);
xt{1} = prior.mu{gs_sel} + chol(prior.sig{gs_sel})' * randn(model.fn,1); 
y_meas{1} = feval(model.hx,xt{1}) + chol(model.R)' * randn(model.hn,1);

% get the trajectory of the sample over time
for k = 2 : time.nSteps
    Ttmp = [time.tspan(k-1) time.tspan(k)];         % integration period   
    xt{k} = sde_int(model, time, Ttmp, xt{k-1});
    y_meas{k} = feval(model.hx,xt{k}) + chol(model.R)' * randn(model.hn,1);
end

%% ------------------------------------------------------------------------
% PF - for particle filters
%--------------------------------------------------------------------------
fprintf('  - create initial particles\n');

% generate i.c. samples
u_pf = rand(1,pf.no_particles);
tu_pf = 0;
X_pf = [];
for i = 1 : prior.n
    % for each gaussian component draw samples dictated by its weight
    % magnitude
    ind = find((u_pf >= tu_pf) & (u_pf < tu_pf + prior.weig(i)));
    Xtmp_pf = repmat(prior.mu{i},1,length(ind)) + chol(prior.sig{i})' * randn(model.fn,length(ind));
    
    % collect all the samples and advance
    X_pf = [X_pf Xtmp_pf];
    tu_pf = tu_pf + prior.weig(i);
end
w_pf = ones(1, pf.no_particles) / pf.no_particles;

%% ------------------------------------------------------------------------
% compute initial estimates
%--------------------------------------------------------------------------
F1.w{1} = prior.weig; 
F2.w{1} = prior.weig; 
F4.w{1} = prior.weig;

[F1.x{1},F1.P{1},F1.ll{1}] =  getGSdata(F1.w{1},prior.mu,prior.sig,X_pf,w_pf);
[F2.x{1},F2.P{1},F2.ll{1}] =  getGSdata(F2.w{1},prior.mu,prior.sig,X_pf,w_pf);
[F4.x{1},F4.P{1},F4.ll{1}] =  getGSdata(F4.w{1},prior.mu,prior.sig,X_pf,w_pf);
[PF.x{1},PF.P{1}] = getPFdata(X_pf,w_pf);

save_weights = prior.weig;

%% ------------------------------------------------------------------------
% init filters
%--------------------------------------------------------------------------
mu = prior.mu;
sig = prior.sig;

save_mus{1}.mus = mu;

%% ------------------------------------------------------------------------
% RUN FILTERS
%--------------------------------------------------------------------------   
for k = 2 : time.nSteps
    fprintf('  - run filters tStep %d / %d\n', k, time.nSteps);    

    % save data 
    mu_old = mu;
    sig_old = sig;
 
    % get the integration interval
    Ttmp = [time.tspan(k-1) time.tspan(k)];
    
%% ------------------------------------------------------------------------
% PF - propagate particles
%--------------------------------------------------------------------------
    fprintf('     - Particle Filter\n'); 
    
    % PF - time update
    X_pf = sde_int(model, time, Ttmp, X_pf);
    
    % PF - measurement update    
    if (measAvail(k) == 1)
        [X_pf, w_pf] = PF_meas_update(X_pf, w_pf, model, pf, y_meas{k});
    end; 
                
    % PF - compute estimates
    [PF.x{k},PF.P{k}] = getPFdata(X_pf, w_pf);
     
%% ------------------------------------------------------------------------
% propagate and update each gaussian component using Kalman Filter
%--------------------------------------------------------------------------
    fprintf('     - Kalman Filter\n');
    
    for j = 1 : prior.n          
        
        % time update
        [mu{j}, sig{j}] = feval([myfilter.name '_time_update'], myfilter, model, mu{j}, sig{j}, Ttmp, opt);
        
        % do update
        if (measAvail(k) == 1)               
            
            [mu_up{j}, sig_up{j}, z_mu{j}, z_sig{j}] = feval([myfilter.name '_measurement_update'], myfilter, model, mu{j}, sig{j}, Ttmp, opt, y_meas{k});         

        else
            mu_up{j} = mu{j};
            sig_up{j} = sig{j};                
        end             
        
    end   
    
%% ------------------------------------------------------------------------
% GS1 - no weight update
%--------------------------------------------------------------------------
    fprintf('     - Weights GS1\n');
    F1.w{k} = F1.w{k-1};    

%% ------------------------------------------------------------------------
% GS2 - Alspach - classic weight update
%--------------------------------------------------------------------------
    fprintf('     - Weights GS2\n');
    if (measAvail(k) == 1)
        for j = 1 : prior.n
            F2.w{k}(j) = F2.w{k-1}(j)*getLikelihood(y_meas{k} - z_mu{j}, z_sig{j} + model.R);
        end
        if (sum(F2.w{k}) > 0)
            F2.w{k} = F2.w{k} ./ sum(F2.w{k});
        else
            F2.w{k} = F2.w{k-1};
        end;
    else
        F2.w{k} = F2.w{k-1}; 
    end 

%% ------------------------------------------------------------------------
% GS4 - Continuous Update - using FPKE
%--------------------------------------------------------------------------   
    fprintf('     - Weights GS4\n');

    switch myfilter.integration_method
        case 'IUT'
            % compute matrix H - Unscented Transformation
            [M,N] = ComputeFPE_QP_UT_noPG(F4.w{k-1}, mu, sig, model, time, k);
        case 'GQ_perComp'
            % compute matrix H - using Gaussian Quadrature with quadrature
            % points for each gaussian component
%             LL = ComputeFPE_QP_GQ_perComp(F4.w{k-1}, mu, sig, model, time, k, myfilter);
%             H = LL+eye(size(LL));            
        case 'GQ_all'
            % compute matrix H - using Gaussian Quadrature with
            % quadrature points that cover the entire domain
%             LL = ComputeFPE_QP_GQ_all(F4.w{k-1}, mu, sig, model, time, k, myfilter);
%             H = LL+eye(size(LL));  
        otherwise
            error('Integration method unknown');
    end;   

    % get the weights
    Aeq = ones(1,prior.n); beq =1; Ain = -eye(prior.n); bin = zeros(prior.n,1);
    xx = sdpvar(prior.n,1);
    xx_old = reshape(F4.w{k-1},prior.n,1);
    cons = set(Aeq*xx == beq) + set(Ain*xx < bin);
	assign(xx, xx_old);
    optyalmip = sdpsettings('solver','sedumi','verbose',0,'usex0',1);
    diagnostic = solvesdp(cons,1/2*xx'*M*xx + xx'*N*xx_old, optyalmip);
    F4.w{k} = double(xx);  
    saveF4w = F4.w{k};

    % do the classic weight update GS2
    if (measAvail(k) == 1)
        for j = 1 : prior.n
            F4.w{k}(j) = F4.w{k}(j)*getLikelihood(y_meas{k} - z_mu{j}, z_sig{j} + model.R);
        end
        if (sum(F4.w{k}) > 0)
            F4.w{k} = F4.w{k} ./ sum(F4.w{k});
        else
            F4.w{k} = saveF4w;
        end;        
    end 
    save_weights = [save_weights; F4.w{k}'];
    F2.w{k}  
    F4.w{k}'
    
%% ------------------------------------------------------------------------
% Store estimates
%--------------------------------------------------------------------------
    mu = mu_up;
    sig = sig_up; 
    
    save_mus{k}.mus = mu;

%% ------------------------------------------------------------------------
% Compute estimates
%--------------------------------------------------------------------------
    [F1.x{k},F1.P{k},F1.ll{k}] =  getGSdata(F1.w{k},mu,sig,X_pf,w_pf);
    [F2.x{k},F2.P{k},F2.ll{k}] =  getGSdata(F2.w{k},mu,sig,X_pf,w_pf);
    [F4.x{k},F4.P{k},F4.ll{k}] =  getGSdata(F4.w{k},mu,sig,X_pf,w_pf);        
   
%% ------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------
    
    plot3(X_pf(1,:),X_pf(2,:),X_pf(3,:),'.r'); hold on;
    for j = 1 : prior.n   
        myX = zeros(model.fn,k);
        for i = 1 : k
            myX(:,i) = save_mus{i}.mus{j};
        end;
        error_ellipse(sig{j},mu{j});
        alpha(F4.w{k}(j));
        plot3(myX(1,:),myX(2,:),myX(3,:),'LineWidth',3);
    end;
    hold off;grid on;view(355,30);
    xlabel('x');ylabel('y');zlabel('z');
    drawnow;
    saveas(gcf,['images/' num2str(k,'%03d')], 'bmp');
    saveas(gcf,['images/' num2str(k,'%03d')], 'fig');    
    pause(0.2);

end