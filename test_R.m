function test_R

L = [ 1;2];
p = [-1;4];

k = 2*(2*rand-1 + 1i*(2*rand-1) );
fprintf('k = %10.8e + i %10.8e\n', real(k), imag(k) );
h = 1e-6*(2*rand-1 + 1i*(2*rand-1) );
fprintf('h = %10.8e + i %10.8e\n', real(h), imag(h) );
N = length(p);

fprintf('Checking make_R and make_dRdk\n');
R = make_R(k,L,p);
C = make_C(k,L,p);
test_C(k,L,p);

z = exp(1i*k*p);
Xns = diag( [1/z(N), z(N)] );
Xs  = [1/2, 1/2i; 1/2, -1/2i];
for n = N:-1:1
    Xn  = diag( [1/z(n), z(n)] );
    Xns = blkdiag(Xn,Xns);
    X   = [1/2, 1/2i; 1/2, -1/2i];
    Xs  = blkdiag(X,Xs);
end

rowscale = eye(2*(N+1));
rowscale(1,1)     = 2*z(1);
rowscale(end,end) = 2/z(N);

abserr = norm(rowscale*C*Xns*Xs - R);
relerr = abserr/norm(R);
fprintf('  R vs C      abserr = %4.2e, relerr = %4.2e\n', abserr, relerr );

Rph = make_R(k+h,L,p);
Rmh = make_R(k-h,L,p);
dRdk = make_dRdk(k,L,p);

abserr = norm( (Rph-Rmh)/(2*h) - dRdk );
relerr = abserr/norm(dRdk);
fprintf('  R vs dRdk   abserr = %4.2e, relerr = %4.2e\n', abserr, relerr );
   
fprintf('Checking make_R_derivs\n');
[Rnew,R_k,R_p,R_kp,R_kk,R_pp] = make_R_derivs(k,L,p);

abserr = norm(R-Rnew);
relerr = abserr/norm(R);
fprintf('  R           abserr = %4.2e, relerr = %4.2e\n', abserr, relerr );
abserr = norm(R_k-dRdk);
relerr = abserr/norm(dRdk);
fprintf('  R_k         abserr = %4.2e, relerr = %4.2e\n', abserr, relerr );

dRdkph = make_dRdk(k+h,L,p);
dRdkmh = make_dRdk(k-h,L,p);

abserr = norm( (dRdkph - dRdkmh)/(2*h) - R_kk );
relerr = abserr/norm(R_kk);
fprintf('  R_kk        abserr = %4.2e, relerr = %4.2e\n', abserr, relerr );

fprintf('  R_p(:,:,n)\n');
for n = 1:length(p)
    en = 0*p; en(n) = 1;
    R1 = make_R(k,L,p+h*en);
    R2 = make_R(k,L,p-h*en);
    R_pn = R_p(:,:,n);
    abserr = norm( (R1-R2)/(2*h) - R_pn );
    relerr = abserr/norm(R_pn);
    fprintf('    n=%d       abserr = %4.2e, relerr = %4.2e\n', ...
                   n,           abserr,         relerr );
end

fprintf('  R_kp(:,:,n)\n');
for n = 1:length(p)
    en = 0*p; en(n) = 1;
    R_k1 = make_dRdk(k,L,p+h*en);
    R_k2 = make_dRdk(k,L,p-h*en);
    R_kpn = R_kp(:,:,n);
    abserr = norm( (R_k1-R_k2)/(2*h) - R_kpn );
    relerr = abserr/norm(R_kpn);
    fprintf('    n=%d       abserr = %4.2e, relerr = %4.2e\n', ...
                   n,           abserr,         relerr );
end

fprintf('  R_pp(:,:,m,n)\n');
for m = 1:length(p)
    for n = 1:length(p)
        en = 0*p; en(n) = 1;
        [~,~,R_p1,~,~,~] = make_R_derivs(k,L,p+h*en);
        [~,~,R_p2,~,~,~] = make_R_derivs(k,L,p-h*en);
        R_p1m = R_p1(:,:,m);
        R_p2m = R_p2(:,:,m);
        R_ppmn = R_pp(:,:,m,n);
        abserr = norm( (R_p1m-R_p2m)/(2*h) - R_ppmn );
        relerr = abserr/norm(R_ppmn);
        fprintf('    m=%d, n=%d  abserr = %4.2e, relerr = %4.2e\n', ...
                       m,    n,       abserr,      relerr );
    end
end

end



function test_C(k,L,p)
% Checks that equations for total/scattered wave set up correctly.

    % --- conditions enforced on total wave -------------------------------

    C = make_C(k,L,p);
    rhs = zeros(length(C),1); rhs(1) = 1;

    % total wave coeffs
    AB = reshape(C\rhs, 2, length(C)/2);
    A = AB(1,:);  B = AB(2,:);
    z = exp(1i*k*p);

    fprintf('  Sanity check compatibility equations (C)\n');
    for ii = 1:length(p)
      left  = A(ii  )*z(ii) + B(ii  )/z(ii);
      right = A(ii+1)*z(ii) + B(ii+1)/z(ii);
      fprintf('    continuity err at delta %d is %g\n', ...
                     ii, abs(right-left) );
      dleft  = 1i*k*(A(ii  )*z(ii) - B(ii  )/z(ii));
      dright = 1i*k*(A(ii+1)*z(ii) - B(ii+1)/z(ii));
      fprintf('    diff       err at delta %d is %g\n', ...
                     ii, abs(dright-dleft - L(ii)*left) );
    end

    fprintf('    A(1)   = %g (should be 1)\n', A(1));
    fprintf('    B(end) = %g (should be 0)\n', B(end));

    % --- C(k) same for scattered and total waves -------------------------

    % compute scattering solution (different rhs)
    rhs = [0*p'; -(L.*z).']; 
    rhs = rhs(:);  
    rhs = [0; rhs; 0];
    ABscatt = reshape(C\rhs, 2, length(C)/2);
    Ascatt = ABscatt(1,:);  Bscatt = ABscatt(2,:);

    fprintf('  Consistency between total and scatt wave comps\n');
    fprintf('    err for ''As'' = %g\n', norm(Ascatt+1 - A));
    fprintf('    err for ''Bs'' = %g\n', norm(Bscatt   - B));

end