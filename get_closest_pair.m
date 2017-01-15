function [u,k] = get_closest_pair(k0,L,p,residtol)
% Finds the eigenpair u,k of R with potential parameters L and p
% with k closest to k0.

if nargin < 4, residtol = 1e-6; end
    
evals = get_resonances(L,p);
evals = filter(evals,L,p,residtol);

% the one closest to k0
[~,ind] = min( abs(evals-k0) );
k = evals(ind);

[u,k] = refine(k,L,p);

end



function [u,k] = refine(k,L,p)
% Newton iteration on bordered system. 
% R(z,L,p) singular iff g(z) zero.

    R = make_R(k,L,p);

    b = rand(length(R),1);
    c = rand(1,length(R));

    for ii = 1:5
        dR = make_dRdk(k,L,p);
        
        tmp = [R, b; c, 0]\[0*b; 1];
        g = tmp(end);

        tmp = -[R, b; c, 0]\[dR*tmp(1:end-1); 0];
        dg = tmp(end);

        dk = -g/dg;
        k = k + dk;
        R  = make_R(k,L,p);
    end
    
    tmp = [R, b; c, 0]\[0*b; 1];
    u = tmp(1:end-1);  u = u/norm(u);
    
end



function evals = filter(evals,L,p,residtol)

    h = p(1:end-1) - p(2:end); maxh = max(abs(h));

    resids = 0*evals;
    for ii = 1:length(evals)
        k  = evals(ii);  absy = abs(imag(k));
      
        if absy > 700/maxh
            resids(ii) = Inf;
        else
            R = make_R(k,L,p);
            resids(ii) = min(svd(R));
        end
    end
    
    evals = evals( resids < residtol );

end



function evals = get_resonances(L,p)

% Pick one:

% --- Eigs from Cheb interp of wave ---------------------------------------
% -------------------------------------------------------------------------

N = 30; % # of interp nodes per interval length 1
evals = resonances_chebsol(L,p,N);


% --- Eigs from Cheb interp of R (in variable k) --------------------------
% -------------------------------------------------------------------------
% To use this one you should already have some idea of a curve
% along which the resonances you're computing are located, and
% then use zspec to map [-1,1] to that curve.

% T = @(z) make_R(z,L,p);
% zspec = [-4-1i; 4-1i];
% N = 50; % # of interp nodes
% evals = cheb_eig(T, zspec, N);


% --- Eigs from Beyn algorithm --------------------------------------------
% -------------------------------------------------------------------------
% This isn't working well as is for this problem.

% % Define boundary of region where we look for the eigs/resonances
% c = -1i;  r = 10;
% phi = @(t) c + r*exp(1i*t);
% dphi = @(t) 1i*r*exp(1i*t);
% inC = @(z) abs(z-c) < r;
% 
% countest = 1;     % Estimate # of resonances in contour phi
% N = 200;          % # points in trapezoid rule eval of integrals
% tolrank = 1e-4;   % below this a singular value considered 0
% tolres  = 1e-1;   % with residual below this eigenvalue accepted
% maxcond = 10;     % above this eigenvalue considered ill-conditioned
% 
% % l and K chosen automatically
% m = length(T(0));
% if countest >= m
%     l = m;
%     K = ceil(countest/m);
% else
%     l = countest;
%     K = 1;
% end
% 
% [evals,~,~,~,~,~,~] ...
%     = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres, ...
%                          maxcond,inC);

end