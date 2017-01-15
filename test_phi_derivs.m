function test_phi_derivs

h = 1e-6;
L = [1;1;2];  
p = [-1;1;2];  

k0 = 12*(rand + 1i*rand);
fprintf('k0 = %10.8e + i %10.8e\n', real(k0), imag(k0) );

residtol = 1e-7;
[u0,k] = get_closest_pair(k0,L,p,residtol);
[hess_phi,grad_phi] = make_phi_derivs(u0,k,k0,L,p);

fprintf('Checking hess_phi and grad_phi errs\n');

N = length(p);

approxgrad_phi = 0*grad_phi;
for n = 1:N
    en = 0*p; en(n) = 1;
    p1 = p + h*en;  phi1 = make_phi(k0,L,p1,residtol);
    p2 = p - h*en;  phi2 = make_phi(k0,L,p2,residtol);
    approxgrad_phi(n) = (phi1 - phi2)/(2*h);
end
abserr = norm(grad_phi - approxgrad_phi);
relerr = abserr/norm(grad_phi);
fprintf('  grad_phi abserr = %4.2e, relerr = %4.2e\n', abserr, relerr );

approxhess_phi = 0*hess_phi;
for m = 1:N
    em = 0*p; em(m) = 1;
    p1 = p + h*em;  [u01,k1] = get_closest_pair(k0,L,p1,residtol);
    [~,grad_phi1] = make_phi_derivs(u01,k1,k0,L,p1);
    p2 = p - h*em;  [u02,k2] = get_closest_pair(k0,L,p2,residtol);
    [~,grad_phi2] = make_phi_derivs(u02,k2,k0,L,p2);
    approxhess_phi(:,m) = (grad_phi1 - grad_phi2)/(2*h);
end
abserr = norm(hess_phi - approxhess_phi);
relerr = abserr/norm(hess_phi);
fprintf('  hess_phi abserr = %4.2e, relerr = %4.2e\n', abserr, relerr );

end



function phi = make_phi(k0,L,p,residtol)

[~,k] = get_closest_pair(k0,L,p,residtol);
phi = abs(k-k0)^2;
    
end
