function C = make_C(k,L,p)
% Computes the matrix of compatibility equations (usually called C).
% Tested within test_R.m. Used for testing R and the matrix in
% resonances_chebsol.m.
%
% Inputs:
%  k - wave number (frequency parameter)
%  L - potential strengths
%  p - potential positions (ordered from left to right)

N  = length(L);
C  = zeros(2*(N+1));

C(1,1)     = 1; % Equation for A(1) (see total_wave)
C(end,end) = 1; % Equation for B(end)

% Compatibility equations at each potential
for n = 1:N
    znk  = exp( 1i*k*p(n));
    znki = exp(-1i*k*p(n));
    Lpik = L(n)/2 + 1i*k;
    Lmik = L(n)/2 - 1i*k;
    
    C(2*n:2*n+1,2*n-1:2*n+2) = ...
        [     znk,       znki,      -znk,      -znki; ...
         Lpik*znk,  Lmik*znki,  Lmik*znk,  Lpik*znki];    
end