% evals = cheb_eig(Tfun, zspec, N)
%
% Find eigenvalues of Tfun(z) near zspec([-1,1]) by solving a polynomial eigenvalue
% problem based on interpolation.  Uses a colleague matrix formulation.
%
% See "Chebyshev interpolation for nonlinear eigenvalue problems" by
% Effenberger and Kressner for how to linearize the Chebyshev
% approximation.
%
% See section 5 of "An extension of MATLAB to continuous functions and 
% operators" by Battles and Trefethen for brief description of how to
% obtain coefficients in cosine series/first-kind Chebyshev expansion (same
% because T_k(cos t) = cos(kt)). Also refers to Spectral Methods in Matlab
% for more details.

function evals = cheb_eig(Tfun, zspec, N)

% -- Default values
if nargin < 2, zspec = [-1, 1]; end
if nargin < 3, N = 10;          end

% -- Form a mapping function from zspec
if isnumeric(zspec)
  zfun = @(x) ( (1-x)*zspec(1) + (1+x)*zspec(2) )/2;
else
  zfun = zspec;
end

% -- Form mapped Chebyshev grid
xn = -cos((0:N)*pi/N).';
zn = zfun(xn);

% -- Set up flattened coefficient array
Tj = Tfun(zn(1));
n  = length(Tj);
Tc = zeros(n*n, N+1);
Tc(:,1) = reshape(Tj, n*n, 1);
for j = 2:N+1
  Tc(:,j) = reshape( Tfun(zn(j)), n*n, 1 );
end

% -- Find coefficients in Chebyshev/cosine expansion
% Note the ordering in input to ifft is due to reordering of Chebyshev
% points above (from left to right rather than usual right to left)
Ac0 = ifft( [fliplr(Tc), Tc(:,2:end-1)], [], 2 );  % DCT would be faster
Ac  = [Ac0(:,1), 2*Ac0(:,2:N), Ac0(:,N+1)];

% -- Set up colleague matrix
e = ones(N-1,1);
C0 = ( diag(e,-1) + diag(e,1) )/2;
C0(1,2) = 1;
C1 = kron(C0, eye(n));
AN = reshape(Ac(:,end), n, n);
Ar = reshape(Ac(:,1:end-1), n, n*N);

C1(end-n+1:end,:) = AN*C1(end-n+1:end,:) - Ar/2;
C2 = eye(length(C1));
C2(end-n+1:end,end-n+1:end) = AN;

% -- Compute eigenvalues and map back to z space
l0 = eig(C1,C2);
evals = zfun(l0);

evals = evals( abs(evals) < Inf );