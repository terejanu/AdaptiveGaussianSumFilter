%--------------------------------------------------------------------------
% Simulations for the Adaptive Weights Gaussian Sum Filter paper
%
% Gabriel Terejanu (terejanu@buffalo.edu)
%--------------------------------------------------------------------------
close all;
clear all;
clc;

%% ------------------------------------------------------------------------
% init
%--------------------------------------------------------------------------
runLocal = true;
mcRuns = 1;

% if running in parallel make sure this are set right
cwd = '/home/csgrad/terejanu/gabi/gs-da/simulations_continuous';
inc = '/home/csgrad/terejanu/matlab_util';

%% ------------------------------------------------------------------------
% model
%--------------------------------------------------------------------------

% ---- LORENZ EXAMPLE -----
model.fn = 3;               % state space dimensionality
model.fx = 'fx_lorenz';
model.hn = 1;               % measurement dimensionality
model.hx = 'hx_lorenz';
% model.Q = diag([0 0 1]); %model.Q(3,3) = 2;
model.Q = diag([1 1 1]); model.Q(3,3) = 2;
model.sQ = sqrt(model.Q);
model.R = 1;

% ---- 2D EXAMPLE -----
% model.fn = 2;               % state space dimensionality
% model.fx = 'fx_2dex';
% model.hn = 1;               % measurement dimensionality
% model.hx = 'hx_2dex';
% model.Q = 0.2^2*ones(2);
% model.sQ = sqrt(model.Q);
% model.R = 1;

% ---- SIN(X) EXAMPLE -----
% model.fn = 1;               % state space dimensionality
% model.fx = 'fx_sinx';
% model.hn = 1;               % measurement dimensionality
% model.hx = 'hx_squared';
% model.Q = 1;
% model.sQ = sqrt(model.Q);
% model.R = 1;

%% ------------------------------------------------------------------------
% prior Gaussian Sum   !!!!!! make them cells (mu & sig) !!!!!!!!
%--------------------------------------------------------------------------

% ---- LORENZ EXAMPLE -----
% prior.n = 2;
% [prior.weig(1), prior.mu{1}, prior.sig{1}, prior.weig(2), prior.mu{2}, prior.sig{2}] = ...
%     split(1, [.5;-1.5;25], eye(3), 1, 0.5);
% prior.weig(1) = 0.5;
% prior.weig(1) = 0.5;

prior.mu{1} = [-.2 -.2 8]';
prior.mu{2} = [.2 .2 8]';
% prior.mu{3} = [0 0 8]';
prior.sig{1} = sqrt(0.35)*eye(3,3);
prior.sig{2} = sqrt(0.35)*eye(3,3);
% prior.sig{3} = sqrt(0.35)*eye(3,3);
prior.n = 2;
prior.weig = [0.9 0.1];
% prior.weig = [0.2 0.2 0.6];

% [X,Y] = meshgrid(1:0.2:2,-3:0.2:-2);
% prior.n = prod(size(X));
% ix = 0;
% for i = 1 : size(X,1)
%     for j = 1 : size(X,2)
%        ix = ix + 1;
%        prior.mu{ix} = [X(i,j) Y(i,j) 6]';
%        prior.sig{ix} = sqrt(0.35)*eye(3,3);
%     end;
% end;
% prior.weig = ones(1,prior.n)/prior.n;

% ---- 2D EXAMPLE -----
% prior.mu{1} = [0.1 0]';
% prior.mu{2} = [-0.1 0]';
% prior.sig{1} = 0.1^2*eye(2);
% prior.sig{2} = 0.1^2*eye(2);
% prior.n = 2;
% prior.weig = [0.9 0.1];

% ---- SIN(X) EXAMPLE -----
% prior.n = 2;
% prior.weig = [0.1 0.9];
% prior.mu{1} = -0.2;
% prior.mu{2} = 0.2;
% prior.sig{1} = 1;
% prior.sig{2} = 1;

%% ------------------------------------------------------------------------
% time
%--------------------------------------------------------------------------

% ---- LORENZ EXAMPLE -----
time.t0 = 0;
time.dt = 0.01;
time.tf = 0.51;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);
time.dt_N = 1000;
time.dt_sde = time.dt / time.dt_N;

% ---- 2D EXAMPLE -----
% time.t0 = 0;
% time.dt = 0.25;
% time.tf = 10;
% time.tspan = time.t0 : time.dt : time.tf;
% time.nSteps = length(time.tspan);

% ---- SIN(X) EXAMPLE -----
% time.t0 = 0;
% time.dt = 0.1;
% time.tf = 8;
% time.tspan = time.t0 : time.dt : time.tf;
% time.nSteps = length(time.tspan);

%% ------------------------------------------------------------------------
% measurements available every n time steps
%--------------------------------------------------------------------------
measTs = 1;
measAvail = zeros(1,time.nSteps);
measAvail(1:measTs:time.nSteps) = 1;

%% ------------------------------------------------------------------------
% particle filter settings
%--------------------------------------------------------------------------
pf.no_particles = 10000;
pf.no_bins = 100;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;                  % treshold for the number of effective samples

pf.plot_pdf = true;             % if mcRuns = 1 makes sense to plot the pdf in time
noPdfs = 10;                     % number pdf slices to be plotted
pf.plot_pdf_available = round(linspace(1,time.nSteps,noPdfs));

pf.bounds = 0.05;               % enlarge the bounds when plotting the pdf (bounds are given by particle filter)

%% ------------------------------------------------------------------------
% filter used for propagation of the 2 momemnts
%--------------------------------------------------------------------------

%----------------------------------------
% avail: EKF    - Extended Kalman Filter
myfilter.name = 'EKF';                  

%----------------------------------------
% avail: IUT            - Unscented Transformation 6n+1 sigma points
myfilter.integration_method = 'IUT';    

%        GQ_perComp     - Gaussian Quadrature - quadrature points obtained
%                         for each gaussian component
% myfilter.integration_method = 'GQ_perComp';  
% myfilter.integration_GQ_bound = 3;  % max bound given the standard deviation of a gaussian component
% myfilter.integration_GQ_no_points = 100; % number of points per dimension

%        GQ_all         - Gaussian Quadrature - quadrature points obtained
%                         for the entire domain
% myfilter.integration_method = 'GQ_all';  
% myfilter.integration_GQ_bound = 3;  % max bound given the standard deviation of a gaussian component
% myfilter.integration_GQ_no_points = 100; % number of points per dimension
                                               
%% ------------------------------------------------------------------------
% init Job Distribution
%--------------------------------------------------------------------------
if (~runLocal)
    jm = findResource('jobmanager','name','underground');

    job = createJob(jm); 
    set(job,'FileDependencies',{cwd}) 

    for iRun = 1 : mcRuns
        obj(iRun) = createTask(job, @runFilters, 5, {iRun,model,prior,time,measAvail,pf,inc,myfilter});
    end;
    submit(job);
    waitForState(job,'finished');
    result = getAllOutputArguments(job);
end;

%% ------------------------------------------------------------------------
% init for MC results
%--------------------------------------------------------------------------

% error wrt the truth
abserrMC_GS1 = zeros(model.fn,length(time.tspan));
abserrMC_GS2 = zeros(model.fn,length(time.tspan));
abserrMC_GS4 = zeros(model.fn,length(time.tspan));
abserrMC_PF = zeros(model.fn,length(time.tspan));

% error wrt the truth
errMC_GS1 = zeros(size(time.tspan));
errMC_GS2 = zeros(size(time.tspan));
errMC_GS4 = zeros(size(time.tspan));
errMC_PF = zeros(size(time.tspan));

% error wrt the particle filter
errPF_GS1 = zeros(size(time.tspan));
errPF_GS2 = zeros(size(time.tspan));
errPF_GS4 = zeros(size(time.tspan));

% average loglikelihood
ll_GS1 = zeros(size(time.tspan));
ll_GS2 = zeros(size(time.tspan));
ll_GS4 = zeros(size(time.tspan));

% average estimator variance / dimension
s_GS1 = zeros(size(time.tspan));
s_GS2 = zeros(size(time.tspan));
s_GS4 = zeros(size(time.tspan));



for iRun = 1 : mcRuns
    
%--------------------------------------------------------------------------
% distributed collection
%--------------------------------------------------------------------------    
    if (~runLocal)
        taskoutput = get(obj(iRun), 'OutputArguments');
        xt = taskoutput{1};     % truth
        F1 = taskoutput{2};     % no update
        F2 = taskoutput{3};     % classic Bayes update
        F4 = taskoutput{4};     % continuous update - FPKE
        PF = taskoutput{5};     % particle filter
    end;
    
%--------------------------------------------------------------------------
% local run
%--------------------------------------------------------------------------       
    if (runLocal)
        [xt,F1,F2,F4,PF] = runFilters(iRun,model,prior,time,measAvail,pf,inc,myfilter);
    end;
    saveALL.xt{iRun} = xt;
    saveALL.F1{iRun} = F1;
    saveALL.F2{iRun} = F2;
    saveALL.F4{iRun} = F4;
    saveALL.PF{iRun} = PF;
    
%--------------------------------------------------------------------------
% get the estimates (1st & 2nd moment)
%--------------------------------------------------------------------------   
    xGS1 = F1.x; xGS2 = F2.x; xGS4 = F4.x; xPF = PF.x;
    sGS1 = F1.P; sGS2 = F2.P;sGS4 = F4.P; sPF = PF.P;
    lGS1 = F1.ll; lGS2 = F2.ll; lGS4 = F4.ll;
    
%--------------------------------------------------------------------------
% get the estimates (1st & 2nd moment, loglikelihood)
%--------------------------------------------------------------------------

%     for i = 1 : model.fn
%         for j = 1 : time.nSteps
%             abserrMC_GS1(i,j) = abserrMC_GS1(i,j) + sqrt(mean((xGS1{j}(i) - xt{j}(i)).^2)) / mcRuns;
%             abserrMC_GS2(i,j) = abserrMC_GS2(i,j) + sqrt(mean((xGS2{j}(i) - xt{j}(i)).^2)) / mcRuns;
%             abserrMC_GS4(i,j) = abserrMC_GS4(i,j) + sqrt(mean((xGS4{j}(i) - xt{j}(i)).^2)) / mcRuns;  
%             abserrMC_PF(i,j) = abserrMC_PF(i,j) + sqrt(mean((xPF{j}(i) - xt{j}(i)).^2)) / mcRuns;            
%         end;
%     end;

    for i = 1 : time.nSteps
        errMC_GS1(i) = errMC_GS1(i) + norm(xGS1{i} - xt{i})^2 / mcRuns;
        errMC_GS2(i) = errMC_GS2(i) + norm(xGS2{i} - xt{i})^2 / mcRuns;
        errMC_GS4(i) = errMC_GS4(i) + norm(xGS4{i} - xt{i})^2 / mcRuns;  
        errMC_PF(i) = errMC_PF(i) + norm(xPF{i} - xt{i})^2 / mcRuns;
        
        errPF_GS1(i) = errPF_GS1(i) + norm(xGS1{i} - xPF{i})^2 / mcRuns;
        errPF_GS2(i) = errPF_GS2(i) + norm(xGS2{i} - xPF{i})^2 / mcRuns;
        errPF_GS4(i) = errPF_GS4(i) + norm(xGS4{i} - xPF{i})^2 / mcRuns;  
        
        ll_GS1(i) = ll_GS1(i) + lGS1{i} / mcRuns;
        ll_GS2(i) = ll_GS2(i) + lGS2{i} / mcRuns;
        ll_GS4(i) = ll_GS4(i) + lGS4{i} / mcRuns;  
        
        s_GS1(i) = s_GS1(i) + trace(sGS1{i}) / model.fn / mcRuns;
        s_GS2(i) = s_GS2(i) + trace(sGS2{i}) / model.fn / mcRuns; 
        s_GS4(i) = s_GS4(i) + trace(sGS4{i}) / model.fn / mcRuns;        
    end;
end        
        
%% ------------------------------------------------------------------------
% destroy job
%--------------------------------------------------------------------------
if (~runLocal)
    destroy(job);
end;

%% ------------------------------------------------------------------------
% PLOT RESULTS
%--------------------------------------------------------------------------

%% - ABS ERR (wrt truth) -----------------------------------------------------
% for i = 1 : model.fn
%     figure; hold on;
%     plot(time.tspan,abserrMC_GS1(i,:),'g');
%     plot(time.tspan,abserrMC_GS2(i,:),'b');
%     plot(time.tspan,abserrMC_GS4(i,:),'r');
%     plot(time.tspan,abserrMC_PF(i,:),'k');
%     legend('GS1 - no update','GS2 - classic','GS4 - w.upd.FPKE','PF - partile filter'); 
%     title(['avg. ABS ERR comparison (wrt truth) - state ' num2str(i)]);
%     hold off;
% end;

%% - RMSE (wrt truth) -----------------------------------------------------
figure; hold on;
plot(time.tspan,sqrt(errMC_GS1),'g');
plot(time.tspan,sqrt(errMC_GS2),'b');
plot(time.tspan,sqrt(errMC_GS4),'r');
plot(time.tspan,sqrt(errMC_PF),'k');
legend('GS1 - no update','GS2 - classic','GS4 - w.upd.FPKE','PF - partile filter'); 
title(['avg. RMSE comparison (wrt truth) - ' num2str(mcRuns)]);
hold off;

%% - RMSE (wrt particle filter) -------------------------------------------
figure; hold on;
plot(time.tspan,sqrt(errPF_GS1),'g');
plot(time.tspan,sqrt(errPF_GS2),'b');
plot(time.tspan,sqrt(errPF_GS4),'r');
legend('GS1 - no update','GS2 - classic','GS4 - w.upd.FPKE'); 
title(['avg. RMSE comparison (wrt particle filter) - ' num2str(mcRuns)]);
hold off;

%% - LOG-LIKELIHOOD -------------------------------------------------------
figure; hold on;
plot(time.tspan,ll_GS1,'g');
plot(time.tspan,ll_GS2,'b');
plot(time.tspan,ll_GS4,'r');
legend('GS1 - no update','GS2 - classic','GS4 - w.upd.FPKE'); 
title(['avg. Loglikelihood of the particles - ' num2str(mcRuns)]);
hold off;

%% - LOG-LIKELIHOOD NORMALIZED --------------------------------------------
% m_div = sqrt(ll_GS2.^2 + ll_GS4.^2);
% 
% figure; hold on;
% plot(time.tspan,ll_GS2./m_div,'b');
% plot(time.tspan,ll_GS4./m_div,'r');
% legend('Classic GSF','Adaptive GSF'); 
% title(['avg. Loglikelihood of the particles - ' num2str(mcRuns)]);
% hold off;

%% - AVERAGE NORMALIZED CROSS-ENTROPY -------------------------------------

nH2 = zeros(size(ll_GS2));
nH4 = zeros(size(ll_GS4));

for i = 1 : size(ll_GS2,1)
    tmp2 = -ll_GS2(i,:)./pf.no_particles;
    tmp4 = -ll_GS4(i,:)./pf.no_particles;
    tmp_div = sqrt(tmp2.^2+tmp4.^2);
       
    if tmp_div > 0,
        nH2(i,:) =  tmp2./tmp_div;
        nH4(i,:) =  tmp4./tmp_div;
    end;
end;

figure; hold on;
plot(time.tspan,mean(nH2,1),'b');
plot(time.tspan,mean(nH4,1),'r');
legend('GS2 - classic','GS4 - w.upd.CKE'); 
title(['Average Normalized Cross-Entropy - ' num2str(mcRuns)]);
hold off;

%% - VARIANCE / DIMENSION -------------------------------------------------
% figure; hold on;
% plot(time.tspan,s_GS1,'g');
% plot(time.tspan,s_GS2,'m');
% plot(time.tspan,s_GS4,'b');
% legend('GS1 - no update','GS2 - classic','GS3 - w.upd.CKE','GS4 - w.upd.FPKE'); 
% title(['avg. Variance per dimension - ' num2str(mcRuns)]);
% hold off;

%% - PLOT PDF SNAPSHOTS ---------------------------------------------------

% if (pf.plot_pdf)
%     
%     % loop over state dimensions
%     for i = 1 : model.fn
%         figure;
%         for k = 1 : length(pf.plot_pdf_available)
%             idx = pf.plot_pdf_available(k);
% 
%             subplot(221); hold on; % - Particle filter
%             patch(PF.pdf{idx}.dim(i).mmX,PF.pdf{idx}.dim(i).mmY,PF.pdf{idx}.dim(i).mmZ,'b');
%             legend('PF');
% 
%             subplot(222); hold on; % - GS2 - classic update
%             patch(F2.pdf{idx}.dim(i).mmX,F2.pdf{idx}.dim(i).mmY,F2.pdf{idx}.dim(i).mmZ,'b');
%             legend('GS2 - classic');
% 
%             subplot(224); hold on; % - GS4 - continuous update - FPKE
%             patch(F4.pdf{idx}.dim(i).mmX,F4.pdf{idx}.dim(i).mmY,F4.pdf{idx}.dim(i).mmZ,'b');
%             legend('GS4 - FPKE');
% 
%         end
%         subplot(221); view(-30,37); grid on; hold off;
%         subplot(222); view(-30,37); grid on; hold off;
%         subplot(224); view(-30,37); grid on; hold off;
%         title('Our comparisson');
% 
% %         figure; hold on; % - GS1 - no update
% %         for k = 1 : length(pf.plot_pdf_available)
% %             idx = pf.plot_pdf_available(k);
% %             patch(F1.pdf{idx}.dim(i).mmX,F1.pdf{idx}.dim(i).mmY,F1.pdf{idx}.dim(i).mmZ,'b');
% %         end
% %         view(-30,37); grid on; hold off;
% %         title('Pdf if no update whatsoever');    
%     end;
%     
% end;

save final_results;