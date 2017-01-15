%
% R = make_R(k,L,p)
%
% Computes the matrix R to be used in residual tests.
%
% Inputs:
%  k - wave number (frequency parameter)
%  L - potential strengths
%  p - potential positions (ordered from left to right)
 
function R = make_R(k,L,p)

N = length(L);
R = zeros(2*(N+1));

R( 1 ,     1:2  ) = [1, -1i];
R(end, end-1:end) = [1,  1i];

h = p(1:end-1)-p(2:end);
s = sin(k*h);
c = cos(k*h);

for n = 1:N-1
    
    Ln = L(n)/2;
    sn = s(n);
    cn = c(n);

    R(2*n:2*n+1,2*n-1:2*n+2) = ...
        [  1,  0,           -cn,           -sn; ...
          Ln,  k,  Ln*cn + k*sn,  Ln*sn - k*cn];    
end

Ln = L(N)/2;

R(2*N:2*N+1,2*N-1:2*N+2) = [ 1,  0,  -1,   0; ...
                            Ln,  k,  Ln,  -k];