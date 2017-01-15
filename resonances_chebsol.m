function evals = resonances_chebsol(L,p,N)
% Approx solution by Cheb interpolants to get approximate eigenvalues

    if nargin < 5, range = [-100, 100, -10, 10]; end
    if nargin < 4, restol = 1e-7; end
    if nargin < 3, N = 30; end

    [M0,M1,M2] = make_M(L,p,N);
%     k = rand + 1i*rand;
%     err = check_M(M0,M1,M2,k,L,p,N)

    evals = polyeig(M0,M1,M2);
    evals = evals( abs(evals) < Inf );
   
end



function [M0,M1,M2] = make_M(L,p,N)

    h     = p(2:end)-p(1:end-1);
    Ns    = ceil( h*N ); 
    sumNs = cumsum([0; Ns]);
    
    M0 = zeros(sum(Ns)); M1 = M0; M2 = M0;

    % differential eqtn between potentials/boundary
    for jj = 1:(length(p)-1)
      Nj = Ns(jj);  I = sumNs(jj) + (1:Nj);  id = eye(Nj);
      [D,~] = cheb(Nj-1);  Dj = D*2/h(jj);  Djsqr = Dj^2;
      M0(I(2:end-1),I) = Djsqr(2:end-1,:);
      M2(I(2:end-1),I) = id(2:end-1,:);
    end

    % continuity and jumps across deltas
    for jj = 2:(length(p)-1)
      I1 = sumNs(jj);         I2 = I1 + 1;
      J1 = sumNs(jj-1)+1:I1;  J2 = J1(end) + (1:Ns(jj));

      M0(I1,J1(end)) =  1;  M0(I1,J2(1))   = -1;

      [D,x] = cheb(Ns(jj-1)-1);  Dj = D*2/h(jj-1);
      M0(I2,J1) = Dj(end,:);  M0(I2,J1(end)) = M0(I2,J1(end)) + L(jj);

      [D,x] = cheb(Ns(jj)-1);    Dj = D*2/h(jj);
      M0(I2,J2) = -Dj(1,:);
    end

    % scattering boundary conditions
    [D,x] = cheb(Ns(1)-1);  Dj = D*2/h(1);
    M0(1,1:Ns(1)) = Dj(1,:);   M0(1,1) = M0(1,1) - L(1);
    M1(1,1) = 1i;

    [D,x] = cheb(Ns(end)-1);  Dj = D*2/h(end);
    M0(end,end-Ns(end)+1:end) = Dj(end,:);
    M0(end,end) = M0(end,end) + L(end);
    M1(end,end) = -1i;

end



function [D,x] = cheb(N)
% CHEB  compute D = differentiation matrix, x = Chebyshev grid
  if N==0, D=0; x=1; return, end
  x = cos(pi*(0:N)/N)'; 
  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
  X = repmat(x,1,N+1);
  dX = X-X';                  
  D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
  D  = D - diag(sum(D'));                 % diagonal entries
  
  % reorder the nodes from left to right
  P = flipud(eye(N+1));
  D = P*D*P;
  x = P*x;
end



function err = check_M(M0,M1,M2,k,L,p,N)

    rhs = zeros(length(M0),1);
    Ns = ceil( (p(2:end)-p(1:end-1))*N );
    
    % continuity and jumps across deltas
    for jj = 2:length(p)-1
      I = sum(Ns(1:jj-1)) + 1;
      rhs(I) = -L(jj)*exp(1i*k*p(jj));
    end
    % scattering boundary condition with jump
    rhs(1)   =  L(1  )*exp(1i*k*p(1  ));
    rhs(end) = -L(end)*exp(1i*k*p(end));

    % solve for scattered wave at interp nodes
    M = M0 + k*M1 + k^2*M2;
    scatt = M\rhs;

    % solve for coeffs for true solution
    C = make_C(k,L,p);
    z = exp(1i*k*p);
    rhs = [0*p'; -(L.*z).']; 
    rhs = rhs(:);  
    rhs = [0; rhs; 0];
    ABscatt = reshape(C\rhs, 2, length(C)/2);
    Ascatt = ABscatt(1,:);  Bscatt = ABscatt(2,:);
    Ascatt = Ascatt(2:end-1);  
    Bscatt = Bscatt(2:end-1);

    % check
    err = 0;
    for jj = 1:length(p)-1
      a = p(jj);  b = p(jj+1);
      [D,x] = cheb(Ns(jj)-1);
      xj = (1-x)*a/2 + (1+x)*b/2;
      truesol = Ascatt(jj)*exp(1i*k*xj) + Bscatt(jj)*exp(-1i*k*xj);
      I = sum(Ns(1:jj-1)) + (1:Ns(jj));
      scattj = scatt(I);
      errj    = norm(truesol - scattj);
      err = err + errj^2;
    end
    err = sqrt(err);

end
