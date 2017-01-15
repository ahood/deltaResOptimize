function [hess_phi,grad_phi] = make_phi_derivs(u0,k,k0,L,p)
% Computes hessian and gradient of phi(p) = abs(k(p)-k0)^2 where k(p) is
% resonance of L,p potential nearest to k0.

% test_k_derivs(1e-6,u0,k,k0,L,p);

[hess_k,grad_k] = make_k_derivs(u0,k,L,p);

% NOTE: phi(p) = real(k(p)-k0)^2 + imag(k(p)-k0)^2
grad_phi = 2*( real(k-k0)*real(grad_k)    + imag(k-k0)*imag(grad_k) );
hess_phi = 2*( real(grad_k)*real(grad_k)' + real(k-k0)*real(hess_k) + ...
               imag(grad_k)*imag(grad_k)' + imag(k-k0)*imag(hess_k) );    

end



function test_k_derivs(h,u0,k,k0,L,p)

    [hess_k,grad_k] = make_k_derivs(u0,k,L,p);
    approxhess_k = 0*hess_k;
    N = length(p);
    
    for n = 1:N
        en = 0*p;  en(n) = 1;
        for m = 1:N
            em = 0*p;  em(m) = 1;
            p1 = p + h*(en+em); [~,k1] = get_closest_pair(k0,L,p1,1e-7);
            p2 = p + h*(en-em); [~,k2] = get_closest_pair(k0,L,p2,1e-7);
            p3 = p + h*(em-en); [~,k3] = get_closest_pair(k0,L,p3,1e-7);
            p4 = p - h*(en+em); [~,k4] = get_closest_pair(k0,L,p4,1e-7);
            approxhess_k(n,m) = (k1-k2-k3+k4)/(4*h^2);
        end
    end
    
    abserr = norm(hess_k - approxhess_k);
    relerr = abserr/norm(hess_k);
    fprintf(' hess_k abserr = %4.2e, relerr = %4.2e\n', abserr, relerr);

end



function [hess_k,grad_k] = make_k_derivs(u0,k,L,p)
%{
Let u0,k be a fixed eigenpair of R(p,k(p)). We impose the
normalization u0'*u(p') = 1 to eigenvectors of R(p',k(p')) 
for the following reason:
if so, then the component of u(p') in the direction of u0 
is always 1, and so derivatives of u with respect to the components
of the position vector p' will always be orthogonal to u0.  

Differentiate R(p,k)u(p) = 0 (wrt p) to get derivatives of k and u.

First deriv wrt p(n):  dnR*u + R*dnu = 0
where
  dnR = (dR/dk)*(dk/dp(n)) + dR/dn (dR/dn deriv wrt n-th arg)
      = R_k    *grad_k(n)  + R_p(:,:,n)
  dnu = du/dp(n)
      = jac_u(:,n)     (RECALL: orthogonal to u0)

This becomes the matrix equation

  [R  R_k*u]*X = -R_p(:,:,n)*u,  X = [jac_u(:,n); grad_k(n)].
 
Then impose

  [u0' 0]*X = 0   (from normalization on eigenvectors u).

          
Second deriv wrt p(n) and p(m): dmdnR*u + dnR*dmu + dmR*dnu + R*dmdnu = 0
where

  dmdnR = dm(dR/dk)*(dk/dp(n)) + (dR/dk)*dm(dk/dp(n)) + dm(dR/dn)
        = dmdRdk   *grad_k(n)  + R_k    *hess_k(n,m)  + dmdRdn
          where
            dmdRdk = (d2R/dk2)*(dk/dp(m)) + dR_k/dm
                   = R_kk     *grad_k(m)  + R_kp(:,:,m)
            dmdRdn = R_kp(:,:,n)*grad_k(m) + R_pp(:,:,n,m)
  dmu   = jac_u(:,m)                 (RECALL: orthogonal to u0)
  dmR   = R_k*grad_k(m) + R_p(:,:,m)
  dmdnu = hess_u(n,m,:)              (RECALL: orthogonal to u0)

This becomes the matrix equation
              
  [R  R_k*u]*Y = b

where

  b = -( (dmdRdk*grad_k(n) + dmdRdn)*u + dnR*jac_u(:,m) + dmR*jac_u(:,n) )

  Y = [hess_u(n,m,:); hess_k(n,m)].

Then impose

  [u0' 0]*Y = 0   (from normalization on eigenvectors u).
%}
    N = length(p);
    jac_u  = zeros(length(u0),N);   % jac_u(:,n) is partial u w.r.t. nth position
    grad_k = zeros(N,1);
    hess_k = zeros(N);

    [R,R_k,R_p,R_kp,R_kk,R_pp] = make_R_derivs(k,L,p);

    A = [ R , R_k*u0; 
         u0',   0   ];

    for n = 1:N
        b = -[R_p(:,:,n)*u0; 0];
        X = A\b;
        jac_u(:,n) = X(1:end-1);
        grad_k(n)  = X(end);
    end

    for n = 1:N
        for m = 1:N
            dmdRdk = R_kk*grad_k(m) + R_kp(:,:,m);
            dmdRdn = R_kp(:,:,n)*grad_k(m) + R_pp(:,:,n,m);
            dmR = R_k*grad_k(m) + R_p(:,:,m);
            dnR = R_k*grad_k(n) + R_p(:,:,n);

            b = -( (dmdRdk*grad_k(n) + dmdRdn)*u0 + ...
                   dnR*jac_u(:,m) + ...
                   dmR*jac_u(:,n) );
            b = [b; 0];

            Y = A\b;
            hess_k(n,m) = Y(end);
        end
    end
    
end
