%
% dRdk = make_dRdk(k,L,p)
%
% Computes the derivative of R(k).
%
% Inputs:
%  k - wave number (frequency parameter)
%  L - potential strengths
%  p - potential positions (ordered from left to right)
 
function dRdk = make_dRdk(k,L,p)

N    = length(L);
dRdk = zeros(2*(N+1));

h = p(1:end-1)-p(2:end);
s = sin(k*h);
c = cos(k*h);

for n = 1:N-1
    
    Ln = L(n)/2;
    sn = s(n);    
    cn = c(n);
    dsn =  cn*h(n);
    dcn = -sn*h(n);

    dRdk(2*n:2*n+1,2*n-1:2*n+2) = ...
        [  0,  0,            -dcn     ,            -dsn     ; ...
           0,  1,  Ln*dcn + sn + k*dsn,  Ln*dsn - cn - k*dcn];    
end

dRdk(2*N:2*N+1,2*N-1:2*N+2) = [ 0,  0,  0,  0; ...
                                0,  1,  0, -1];