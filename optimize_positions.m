%
% [alpha,k] = optimize_positions(k0,L,p)
%
% Using Newton's method, we minimize the function 
% phi(p) = (k(p) - k0)^2
% where k(p) is the resonance closest to the target k0 when the delta
% potentials are at p.
%
% Inputs:
%  k0 - complex number that we want resonance close to
%  p  - ordered vector of initial potential positions
%  L  - corresponding potential strengths
%
% Outputs:
%  alpha - final potential positions
%  k     - resonance found closest to k0

function [alpha,k] = optimize_positions(k0,L,p)

% make sure they're columns
L = L(:);
p = p(:);

N      = length(p);
% Left-right shifts don't affect resonances. WLOG
% let's put the left-most at 0.
alpha  = p;
s = ones(size(alpha));

% start Newton to find zero of phi'
count1 = 0;
while s'*s > 1e-6
    count1 = count1 + 1;
    if count1 >= 50
        fprintf('Tried 50 steps of Newton without achieving convergence.\n');
        break
    end
    [u,k] = get_closest_pair(k0,L,alpha);
    [hess_phi,grad_phi] = make_phi_derivs(u,k,k0,L,alpha);

    % Make hess_phi symmetrical
    hess_phi = (hess_phi + hess_phi')/2;

    % Fix first component (DSB)
    grad_phi = grad_phi(2:end);
    hess_phi = hess_phi(2:end,2:end);
    try
        chol(hess_phi);
    catch
        mu       = -min(eig(hess_phi)) + h;
        hess_phi = hess_phi + mu*eye(N-1);
    end
        
    s = [0; -hess_phi\grad_phi];
    
    % Checking if s is a descent direction
    if 1 
        fprintf('s is ');
        if grad_phi'*s(2:end) >= 0
            fprintf('not ');
        end
        fprintf('a descent direction\n');
    end
    
    alpha1 = alpha + s; % Newton step
        
    % Check if full Newton step works and adjust if necessary
    [u1,k1] = get_closest_pair(k0,L,alpha1);
    count2 = 1;
    while (abs(k1-k0)^2 > abs(k-k0)^2)
        alpha1  = (alpha + alpha1)/2;
        [u1,k1] = get_closest_pair(k0,L,alpha1);
        count2  = count2+1;
        if count2 >= 15
            fprintf('Did 15 steps of line search--giving up now.\n');
            break
        end
    end

    alpha = alpha1;
end

[~,k] = get_closest_pair(k0,L,alpha);
