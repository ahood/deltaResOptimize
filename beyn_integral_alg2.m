% function [evals,evecs,Vhat,B0,B1,k,condin] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC)
%
% The values of certain contour integrals are used to make a constant
% matrix whose eigenstructure matches that of T restricted to the interior
% of Jordan curve C parametrized by phi, in the sense that each have the
% same eigenvalues (including multiplicity) and some other conditions.
% This algorithm allows for k >= m, where k is the number of eigenvalues
% of T in the contour C.
% 
% Inputs:
%   l,K - initial index choices
%   T   - holomorphic m x m matrix-valued function
%   phi - parametrization of Jordan curve
%   phi - derivative of phi
%   N   - number of subintervals in trapezoid rule evaluation of Ap's
%   tolrank  - threshold for considering singular values zero
%   tolres   - residual threshold for accepting eigenpair
%   maxcond  - the largest a condition number can be before the eigenvalue
%              is considered ill-conditioned
%   inC - function to determine whether given point is in C
%
% Outputs:
%   evals  - eigenvalues inside contour
%   evecs  - one eigenvector per eigenvalue
%   Vhat   - the random (or identity) matrix used
%   B0,B1  - matrices of computed integrals
%   k      - computed number of eigs in curve C
%   condin - condition numbers of computed eigenvalues in C (may be more
%            than length of evals, since evals has only one rep from each
%            block)

function [evals,evecs,Vhat,B0,B1,k,condin] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC)

evals = [];  evecs = [];  condin = []; 
m = length(T(0));
t = linspace(0,2*pi,N+1);

k = K*l; % so we do calculations at least once.

iter = 0;  maxiter = 10; % preclude infinite loop
while (k == K*l && iter < maxiter)
  iter = iter + 1;
 
  if l < m
    l = l + 1;
    Vhat = rand(m,l);
    [D,k,V0top,B0,B1] = make_D(Vhat,K,N,T,phi,dphi,t,tolrank);
  else
    K = K + 1;
    Vhat = eye(m);
    [D,k,V0top,B0,B1] = make_D(Vhat,K,N,T,phi,dphi,t,tolrank);
  end
  if k == 0
    return;
  end
  
%   If k = K*l
%      start over and increase l or K by 1. In the remarks, it is said
%      that this is used as an indication that there might be more than K*l 
%      eigs encircled by C.
end

[S,E,c] = condeig(D);

% If all eigs well conditioned
if max(c) < maxcond
  for ii = 1:length(E)
    e = E(ii,ii);  s = S(:,ii);  v = V0top*s;
    if (inC(e) && norm(T(e)*v) < tolres)
      evals = [evals, e];  evecs = [evecs, v/norm(v)];
    end
  end
else
  % Compute Schur decomp B*Q = Q*U ordering U so that eigenvalues
  %   actually encircled by C appear first and discard the rest and
  %   corresponding columns of Q.
  [Q,U] = schur(D,'complex');
  [clusters,BLKS,condin] = make_clusters(U,inC,maxcond);
  [Q,U] = ordschur(Q,U,clusters); % now outside eigs trailing and rest grouped
  n = sum(BLKS); % # eigs in C
  Qin = Q(:,1:n);  Uin = U(1:n,1:n);
  % Block diagonalize U such that diagonal blocks belong to different 
  %   eigenvalues. 
  [Qblk,Ublk] = bdschur(Uin,[],BLKS);  
  % Let [e] be the diagonal entry of [a given] block and determine the 
  %   corresponding eigenvector [s] of B from the first column of the 
  %   [same] block in U. (Evec of Ublk is e_start where start is the column
  %   at which the current block starts.)
  e = Ublk(1,1);
  s = Qin*Qblk(:,1);  v = V0top*s;
  if norm(T(e)*v) < tolres
    evals = [evals, e];  evecs = [evecs, v/norm(v)];
  end
  start = 1;
  for ii = 2:length(BLKS)
    start = start + BLKS(ii-1); % index at which ii'th block starts
    e = Ublk(start,start);
    s = Qin*Qblk(:,start);  v = V0top*s;
    if norm(T(e)*v) < tolres
      evals = [evals, e];  evecs = [evecs, v/norm(v)];
    end
  end
end

end

% Make the matrix B which has the same eigenstructure as T. B and V0 are 
% only computed if k < l. Other outputs are used for testing.

function [D,k,V0top,B0,B1] = make_D(Vhat,K,N,T,phi,dphi,t,tolrank)
  k = 0; D = 0; V0top = 0;
  [m,l] = size(Vhat);

  Ap = zeros(m,l*2*K); % [A0, A1, A2, ..., A_(2*K-1)]
  for ii = 1:N
    temp = T(phi(t(ii)))\Vhat*dphi(t(ii));
    for p = 0:2*K-1
      cols = p*l+1:(p+1)*l;
      Ap(:,cols) = Ap(:,cols) + temp*phi(t(ii))^p;
    end
  end
  Ap = Ap/(1i*N);
  
  B0 = []; B1 = [];
  for ii = 0:K-1
    cols0 = ii*l+1:(ii+K)*l;
    cols1 = cols0 + l;
    B0 = [ B0; Ap(:,cols0) ];
    B1 = [ B1; Ap(:,cols1) ];
  end  
    
  [V,Sig,W] = svd(B0,0);
  ind = find(diag(Sig) > tolrank); 
  
  if ~isempty(ind)
    k = ind(end);
    if k < K*l
      V0 = V(1:K*m,1:k);  Sig0 = Sig(1:k,1:k);  W0 = W(1:K*l,1:k);
      D = V0'*B1*W0/Sig0;
      V0top = V0(1:m,:);
    end
  end
  
end


% function [clusters,BLKS,condin] = make_clusters(U,inC,maxcond)
%
% Takes upper triangular U, associates eigenvalues not in C 
% to the lowest cluster number, and organizes the remaining eigenvalues 
% so that badly conditioned eigenvalues that are sufficiently close will 
% be given the same cluster number. 
%
% Inputs:
%   U       - upper triangular matrix
%   inC     - indicator function for interior of curve C
%   maxcond - if eigenvalue condition number > maxcond then we consider
%             that eigenvalue ill-conditioned
%
% Outputs:
%   clusters - vector of cluster numbers 
%   BLKS     - vector of block sizes *only* for eigs inside C.
%   condin   - vector of condition numbers for the returned eigenvalues

function [clusters,BLKS,condin] = make_clusters(U,inC,maxcond)
  u = diag(U);
  clusters = 0*u;
  dealtwith = logical(0*u); % keeping track of which have been assoc to cluster
  [V,E,c] = condeig(U); % assume order preserved for upper tri (by exper.)
  
  % outside eigs
  indout = logical(0*u); % which are outside C
  for ii = 1:length(indout);
      indout(ii) = ~inC(u(ii));
  end
  numout = sum(indout);
  clusters(indout) = ones(numout,1); % outside eigs last
  dealtwith(indout) = true(numout,1); % done with outside eigs
  
  % well-conditioned eigs
  indinwc = c < maxcond & ~dealtwith;
  numinwc = sum(indinwc);
  clusters(indinwc) = 2*ones(numinwc,1); % wc eigs second to last
  BLKS = ones(1, numinwc); % wc eigs have blocks of size 1
  dealtwith(indinwc) = ones(numinwc,1); % done with wc eigs
  
  % badly-conditioned eigs
  currclusternum = 2;
  while sum(~dealtwith) > 0
    currclusternum = currclusternum + 1;
    currlist = u(~dealtwith);
    
    clustertol = 1e-6; % cluster tolerance for matching to currlist(1) - should be based on condition number and machine epsilon???
    ind = abs( currlist(1) - currlist ) < clustertol;
    currclustersize = sum(ind);
    
    temp = find(dealtwith == 0);
    
    clusters(temp(ind)) = currclusternum*ones(currclustersize,1);
    BLKS = [currclustersize, BLKS];
        
    dealtwith(temp(ind)) = ones(currclustersize,1);
  end
  
  % reorder condition numbers the same way and leave out outside eig ones
  condin = [];
  for ii = 2:currclusternum
    indii = find(clusters == ii);
    condin = [c(indii); condin];
  end
      
end