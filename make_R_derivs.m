function [R,R_k,R_p,R_kp,R_kk,R_pp] = make_R_derivs(k,L,p)

N    = length(p);   
    
R    = zeros(2*(N+1));
R_k  = zeros(2*(N+1));
R_kk = zeros(2*(N+1));
R_p  = zeros(2*(N+1),2*(N+1),N);
R_kp = zeros(2*(N+1),2*(N+1),N);
R_pp = zeros(2*(N+1),2*(N+1),N,N);

R( 1 ,     1:2  ) = [1, -1i];
R(end, end-1:end) = [1,  1i];

h = p(1:end-1)-p(2:end);
s = sin(k*h);
c = cos(k*h);

for n = 1:N-1

    rows = 2*n  :2*n+1;
    cols = 2*n-1:2*n+2;

    Ln = L(n)/2;  hn = h(n);
    
    sn    = s(n);     cn    = c(n);     
    dksn  =   cn*hn;  dkcn  =   -sn*hn;
    ddksn = dkcn*hn;  ddkcn = -dksn*hn;

    R(rows,cols) = ...
        [  1,  0,           -cn,           -sn; ...
          Ln,  k,  Ln*cn + k*sn,  Ln*sn - k*cn]; 

    R_k(rows,cols) = ...
        [  0,  0,            -dkcn      ,            -dksn      ; ...
           0,  1,  Ln*dkcn + sn + k*dksn,  Ln*dksn - cn - k*dkcn];    

    R_kk(rows,cols(3:4)) = ...
        [             -ddkcn         ,             -ddksn         ; ...
          Ln*ddkcn + 2*dksn + k*ddksn, Ln*ddksn - 2*dkcn - k*ddkcn];      
      
    dph = 1; % first derivative wrt p(n)  
    dpsn  =   cn*k*dph;  dpcn  =   -sn*k*dph;
    dpdksn =  dpcn*hn + cn*dph;
    dpdkcn = -dpsn*hn - sn*dph;
      
    R_p(rows,cols(3:4),n) = ...
        [     -dpcn      ,      -dpsn      ; ...
         Ln*dpcn + k*dpsn, Ln*dpsn - k*dpcn];

    R_kp(rows,cols(3:4),n) = ...
        [     -dpdkcn               ,      -dpdksn               ; ...
         Ln*dpdkcn + dpsn + k*dpdksn, Ln*dpdksn - dpcn - k*dpdkcn];

        dph = 1; % second derivative wrt p(n) 
        ddpsn = dpcn*k*dph;  ddpcn = -dpsn*k*dph;

        R_pp(rows,cols(3:4),n,n) = ...
            [     -ddpcn       ,      -ddpsn       ; ...
             Ln*ddpcn + k*ddpsn, Ln*ddpsn - k*ddpcn];

        dph = -1; % second derivative wrt p(n+1) 
        ddpsn = dpcn*k*dph;  ddpcn = -dpsn*k*dph;

        R_pp(rows,cols(3:4),n,n+1) = ...
            [     -ddpcn       ,      -ddpsn       ; ...
             Ln*ddpcn + k*ddpsn, Ln*ddpsn - k*ddpcn];
    
     
    dph = -1; % first derivative wrt p(n+1)
    dpsn  =   cn*k*dph;  dpcn  =   -sn*k*dph;
    dpdksn =  dpcn*hn + cn*dph;
    dpdkcn = -dpsn*hn - sn*dph;
    
    R_p(rows,cols(3:4),n+1) = ...
        [     -dpcn      ,      -dpsn      ; ...
         Ln*dpcn + k*dpsn, Ln*dpsn - k*dpcn];

    R_kp(rows,cols(3:4),n+1) = ...
        [     -dpdkcn               ,      -dpdksn               ; ...
         Ln*dpdkcn + dpsn + k*dpdksn, Ln*dpdksn - dpcn - k*dpdkcn];

        dph = 1; % second derivative wrt p(n) 
        ddpsn = dpcn*k*dph;  ddpcn = -dpsn*k*dph;

        R_pp(rows,cols(3:4),n+1,n) = ...
            [     -ddpcn       ,      -ddpsn       ; ...
             Ln*ddpcn + k*ddpsn, Ln*ddpsn - k*ddpcn];

        dph = -1; % second derivative wrt p(n+1) 
        ddpsn = dpcn*k*dph;  ddpcn = -dpsn*k*dph;

        R_pp(rows,cols(3:4),n+1,n+1) = ...
            [     -ddpcn       ,      -ddpsn       ; ...
             Ln*ddpcn + k*ddpsn, Ln*ddpsn - k*ddpcn];
end

Ln = L(N)/2;

R(2*N:2*N+1,2*N-1:2*N+2)   = [ 1,  0,  -1,   0; ...
                              Ln,  k,  Ln,  -k];
                        
R_k(2*N:2*N+1,2*N-1:2*N+2) = [ 0,  0,   0,   0; ...
                               0,  1,   0,  -1];

end