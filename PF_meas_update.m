function [X_new, w_new] = PF_meas_update(X_old, w_old, model, pf, ymeas)
%
% Particle Filter Measurement Update
%
% X_new   - the new set of particles after reweighting and resampling if necessary
%
% Gabriel Terejanu (terejanu@buffalo.edu)

%% ------------------------------------------------------------------------
% weights update
%--------------------------------------------------------------------------
mX_pf = zeros(model.hn, pf.no_particles);
w_new = zeros(size(w_old));
for j = 1 : pf.no_particles
    % get the measurement
    mX_pf(:,j) = feval(model.hx, X_old(:,j));
    w_new(j) = w_old(j) * getLikelihood(ymeas - mX_pf(:,j), model.R);
end
w_new = w_new./sum(w_new);
X_new = X_old;

%% ------------------------------------------------------------------------
% resampling
%--------------------------------------------------------------------------
if (pf.resample) 
	Neff = 1/sum(w_new.^2);
    if (Neff < pf.neff * pf.no_particles)
        I = myresampling(w_new);
        I = round(I);
        for j = 1 : pf.no_particles
            X_new(:,j) = X_old(:,I(j));
        end;    
        % reset weights
        w_new = ones(1,pf.no_particles)/pf.no_particles;
    end;
end;

