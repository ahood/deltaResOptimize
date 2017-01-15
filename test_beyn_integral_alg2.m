function test_beyn_integral_alg2(verbose)

close all
if nargin == 0, verbose = 0; end

% make sure it passes all alg1 tests

%% Test 1: simple eigenvalues
% T(z) = A-zI, A = diag([1,2])

A = diag([1,2]);  
T = @(z) A - z*eye(2);

% Test contour 1 - contour encloses one of two eigenvalues
c = A(2,2);  r = abs(A(1,1)-A(2,2))/2;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 100;
tolrank = 1e-4;
tolres = 1e-4;
countest = 1;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing T(z) = A-zI, B(%g,%g), N = %d\n\n', c, r, N);
fprintf('    A =   1  0\n');
fprintf('          0  2\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs, Vhat, B0, B1, k, ~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
B0err = norm(B0 - [0*Vhat(1,:); -1*Vhat(2,:)]);
B1err = norm(B1 - [0*Vhat(1,:); -2*Vhat(2,:)]);
evalserr = abs(2 - evals);
evecserr = abs([1;0]'*evecs);

print_results(verbose,B0err,B1err,evalserr,evecserr,k);

% Test contour 2 - contour encloses no eigenvalues
c = 5;  r = 1/2;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 100;
tolrank = 1e-4;
tolres = 1e-4;
countest = 1;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing T(z) = A-zI, B(%g,%g), N = %d\n\n', c, r, N);
fprintf('    A =   1  0\n');
fprintf('          0  2\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs, Vhat, B0, B1, k, ~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
B0err = norm(B0 - 0*Vhat);
B1err = norm(B1 - 0*Vhat);
evalserr = norm(evals);
evecserr = norm(evecs);

print_results(verbose,B0err,B1err,evalserr,evecserr,k);

% Test contour 3 - contour very near an eigenvalue
c = A(2,2);  r = abs(A(2,2)-A(1,1))-1e-2;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 2000;
tolrank = 1e-4;
tolres = 1e-4;
countest = 1;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing T(z) = A-zI, B(%g,%g), N = %d\n\n', c, r, N);
fprintf('    A =   1  0\n');
fprintf('          0  2\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs, Vhat, B0, B1, k, ~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
B0err = norm(B0 - [0*Vhat(1,:); -1*Vhat(2,:)]);
B1err = norm(B1 - [0*Vhat(1,:); -2*Vhat(2,:)]);
evalserr = abs(2 - evals);
evecserr = abs([1;0]'*evecs);

print_results(verbose,B0err,B1err,evalserr,evecserr,k);

%% Test 2: Ill-conditioned eigenvalues
% T(z) = A-zI, A = [1, 1, 0; 0, 1, 0; 0, 0, 2]

A = [1, 1, 0; ...
     0, 1, 0; ...
     0, 0, 2];
T = @(z) A - z*eye(3);

% Test contour 1 - contour encloses simple eigenvalue
c = A(3,3);  r = abs(A(2,2)-A(3,3))/2;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 100;
tolrank = 1e-4;
tolres = 1e-4;
countest = 1;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing T(z) = A-zI, B(%g,%g), N = %d\n\n', c, r, N);
fprintf('          1  1  0\n');
fprintf('    A =   0  1  0\n');
fprintf('          0  0  2\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs, Vhat, B0, B1, k, ~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
B0err = norm(B0 - [0*Vhat(1:2,:); -1*Vhat(3,:)]);
B1err = norm(B1 - [0*Vhat(1:2,:); -2*Vhat(3,:)]);
evalerr = abs(2 - evals);
evecserr = norm(abs(evecs)-[0;0;1]);

print_results(verbose,B0err,B1err,evalserr,evecserr,k);

% Test contour 2 - contour encloses double eigenvalue (geo mult < alg mult)

c = A(2,2);  r = abs(A(2,2)-A(3,3))/2;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 100;
tolrank = 1e-4;
tolres = 1e-4;
countest = 2;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing T(z) = A-zI, B(%g,%g), N = %d\n\n', c, r, N);
fprintf('          1  1  0\n');
fprintf('    A =   0  1  0\n');
fprintf('          0  0  2\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs, Vhat, B0, B1, k, ~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
B0err = norm(B0 - [-1*Vhat(1:2,:); 0*Vhat(3,:)]);
B1err = norm(B1 - [-1*Vhat(1,:) + -1*Vhat(2,:); -1*Vhat(2,:); 0*Vhat(3,:)]);
evalerr = abs(1 - evals);
evecserr = norm(abs(evecs)-[1;0;0]);

print_results(verbose,B0err,B1err,evalserr,evecserr,k);

%% Test 3: both well-conditioned and ill-conditioned enclosed
% T(z) = A - zI, A = [1, 1, 0, 0; 0, 1, 0, 0; 0, 0, 1.1, 0; 0, 0, 0, 2]

A = [1, 1, 0,   0; ...
     0, 1, 0,   0; ...
     0, 0, 1.1, 0; ...
     0, 0, 0,   2];
T = @(z) A - z*eye(4);

c = A(2,2);  r = abs(A(2,2)-A(4,4))/2;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 100;
tolrank = 1e-4;
tolres = 1e-4;
countest = 3;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing T(z) = A-zI, B(%g,%g), N = %d\n\n', c, r, N);
fprintf('          1  1    0  0\n');
fprintf('    A =   0  1    0  0\n');
fprintf('          0  0  1.1  0\n');
fprintf('          0  0    0  2\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs, Vhat, B0, B1, k, ~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
B0err = norm(B0 - [-1*Vhat(1:3,:); 0*Vhat(4,:)]);
B1err = norm(B1 - [-1*Vhat(1,:) + -1*Vhat(2,:); -1*Vhat(2,:); -1.1*Vhat(3,:); 0*Vhat(4,:)]);
evalerr = norm(sort(evals)-sort([1,1.1]));
evecserr = norm(abs(evecs)-[1, 0; 0, 0; 0, 1; 0, 0]);

print_results(verbose,B0err,B1err,evalserr,evecserr,k);

% Now triple eigenvalue with only one well-conditioned

A = [1, 1, 0, 0; ...
     0, 1, 0, 0; ...
     0, 0, 1, 0; ...
     0, 0, 0, 2];
T = @(z) A - z*eye(4);

c = A(2,2);  r = abs(A(2,2)-A(4,4))/2;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 100;
tolrank = 1e-4;
tolres = 1e-4;
countest = 3;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing T(z) = A-zI, B(%g,%g), N = %d\n\n', c, r, N);
fprintf('          1  1  0  0\n');
fprintf('    A =   0  1  0  0\n');
fprintf('          0  0  1  0\n');
fprintf('          0  0  0  2\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs, Vhat, B0, B1, k, condin] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
B0err = norm(B0 - [-1*Vhat(1:3,:); 0*Vhat(4,:)]);
B1err = norm(B1 - [-1*Vhat(1,:) + -1*Vhat(2,:); -1*Vhat(2,:); -1*Vhat(3,:); 0*Vhat(4,:)]);
evalerr = norm(sort(evals)-sort([1,1]));
if size(evecs,2) == 2
  evecserr = norm(evecs - [ [evecs(1,1);0;0;0], [evecs(1,2);0;evecs(3,2);0] ]);
  print_results(verbose,B0err,B1err,evalerr,evecserr,k);
elseif min(condin) > maxcond
  fprintf('  ***  Poor choice of maxcond results in too few eigenvectors  ***\n\n');
else
  fprintf('  ***  Too few eigenvectors, reason unknown  ***\n\n');
end

%% Test 4: Example 4.9 from Beyn's paper

T0 = rand(60); T1 = rand(60); T2 = rand(60);
T = @(z) T0 + z*T1 + T2*(z^2);

R = 0.33;
phi = @(t) R*exp(1i*t);
dphi = @(t) 1i*R*exp(1i*t);
N = 150; % from caption
tolrank = 1e-4;
tolres = 1e-1;
countest = 1;
maxcond = 10;
inC = @(z) abs(z) < R;

fprintf('Testing Beyn Example 4.9\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs,~,~,~,~,~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
[X,e] = polyeig(T0,T1,T2);
fprintf('  See figure.\n\n');
figure; hold on
plot(phi(linspace(0,2*pi,100)));
plot(real(e),imag(e),'*k');
plot(real(evals),imag(evals),'ob'); hold off
axis([-R R -R R])
axis square
title('Test 4: Example 4.9')

%% Test 5: Example 4.11 from Beyn's paper

m = 400;  T1 = diag( 2*ones(m  ,1), 0) + ...
               diag(-1*ones(m-1,1),-1) + ...
               diag(-1*ones(m-1,1), 1);  T1(end,end) = 1;  T1 = T1*m;
          T3 = diag( 4*ones(m  ,1), 0) + ...
               diag(   ones(m-1,1),-1) + ...
               diag(   ones(m-1,1), 1);  T3(end,end) = 2;  T3 = T3/(6*m);
          T2 = 0*T1;                     T2(end,end) = 1;
T = @(z) T1 + T2/(1-z) - z*T3;

c = 150; r = 148;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 50;         % residuals smallest around here (Figure 4)
tolrank = 1e-4; % not specified so I used value from Example 4.9
tolres = 1e-1;  % "
countest = 5;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing Beyn Example 4.11\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs,~,~,~,~,~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);

fprintf('  Output of algorithm:\n');
for ii = 1:length(evals)
    res = norm(T(evals(ii))*evecs(:,ii));
    fprintf('    %+8.6e + i %+8.6e, %8.6e\n', real(evals(ii)), imag(evals(ii)), res);
end
fprintf('\n');
if length(evals) ~= 5
    fprintf('  ***  There are supposed to be 5 eigenvalues  ***\n\n');
end

%% Test 6: Example 4.12/5.4 from Beyn's paper (rank deficiency)

T0 = rand(15); T1 = rand(15);
T0(:,1) = 0*T0(:,1);
a = -0.2;  b = 0.1;
T = @(z) T0 + (z-a)*(b-z)*T1;

R = 0.33;
phi = @(t) R*exp(1i*t);
dphi = @(t) 1i*R*exp(1i*t);
N = 150; % accuracy discussed for this N
tolrank = 1e-4;
tolres = 1e-1;
% countest = 4;
maxcond = 10;
inC = @(z) abs(z) < R;

fprintf('Testing Beyn Example 4.12/5.4\n\n');

[X,e] = polyeig(T0-a*b*T1, (a+b)*T1, -T1);

countest = sum( abs(e) < R );
% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 2;  % handles rank deficiency for this problem
end

[evals, evecs,~,~,~,~,~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
e = sort(e);
evals = sort(evals);
fprintf('  Output of polyeig:\n');
for ii = 1:length(e)
    this = e(ii);
    if norm([real(this),imag(this)],1) < R
        fprintf('    %+8.6e + i %+8.6e\n',real(this),imag(this));
    end
end
fprintf('  Output of algorithm:\n');
for ii = 1:length(evals)
    fprintf('    %+8.6e + i %+8.6e\n',real(evals(ii)),imag(evals(ii)));
end
fprintf('\n');
fprintf('  With K = 2 we do not miss %8.6e and %8.6e\n\n',a,b);

% Case where k == m

%% Test 7: simple eigenvalues with k == m
% T(z) = A-zI, A = diag([1,2])

A = diag([1,2]);  
T = @(z) A - z*eye(2);

% contour encloses both eigenvalues
c = A(2,2);  r = abs(A(1,1)-A(2,2))*2;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 100;
tolrank = 1e-4;
tolres = 1e-4;
countest = 2;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing T(z) = A-zI, B(%g,%g), N = %d\n\n', c, r, N);
fprintf('    A =   1  0\n');
fprintf('          0  2\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs, Vhat, B0, B1, k, ~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
[evals,perm] = sort(evals);
evecs = evecs(:,perm);
% B0err = norm(B0 - [0*Vhat(1,:); -1*Vhat(2,:)]);
% B1err = norm(B1 - [0*Vhat(1,:); -2*Vhat(2,:)]);
B0err = 0; B1err = 0;
evalserr = norm([1,2] - evals);
evecserr = norm([evecs(1,1), 0; 0, evecs(2,2)] - evecs);

print_results(verbose,B0err,B1err,evalserr,evecserr,k);

% k > m

%% Test 8: simple eigenvalues with k > m
% T(z) = diag([(1-z)*(1.5-z), (2-z)])

T = @(z) diag([(1-z)*(1.5-z),2-z]);

% contour encloses all three eigenvalues
c = A(2,2);  r = abs(A(1,1)-A(2,2))*2;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 100;
tolrank = 1e-4;
tolres = 1e-4;
countest = 3;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing B(%g,%g), N = %d\n\n', c, r, N);
fprintf('    T(z) =  (1-z)(1.5-z)   0   \n');
fprintf('                 0        2-z\n\n');

% Choose l and K
m = length(T(0));
if countest >= m
    l = m;
    K = ceil(countest/m);
else
    l = countest; % expect k can only be as large as # eigs in contour
    K = 1;
end
[evals, evecs, Vhat, B0, B1, k, ~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
[evals,perm] = sort(evals);
evecs = evecs(:,perm);
% B0err = norm(B0 - [0*Vhat(1,:); -1*Vhat(2,:)]);
% B1err = norm(B1 - [0*Vhat(1,:); -2*Vhat(2,:)]);
B0err = 0; B1err = 0;
evalserr = norm([1,1.5,2] - evals);
evecserr = norm([evecs(1,1), evecs(1,2), 0; 0, 0, evecs(2,3)] - evecs);

print_results(verbose,B0err,B1err,evalserr,evecserr,k);

%% Test 9: Example 5.5 from Beyn's paper (time delay)

T0 = [-5, 1; 2, -6];  T1 = [-2, 1; 4, -1];  tau = 1;
T = @(z) z*eye(2) - T0 - T1*exp(-z*tau);

c = -1;  r = 6;
phi = @(t) c + r*exp(1i*t);
dphi = @(t) 1i*r*exp(1i*t);
N = 150;
tolrank = 1e-4;
tolres = 1e-4;
countest = 2;
maxcond = 10;
inC = @(z) abs(z-c) < r;

fprintf('Testing Beyn Example 5.5\n\n');

% Choose l and K
% m = length(T(0));
% if countest >= m
%     l = m;
%     K = ceil(countest/m);
% else
%     l = countest; % expect k can only be as large as # eigs in contour
%     K = 1;
% end
l = 2; K = 3;
[evals, evecs, Vhat, B0, B1, k, ~] = beyn_integral_alg2(l,K,T,phi,dphi,N,tolrank,tolres,maxcond,inC);
if length(evals) == 5
    fprintf('  5 eigenvalues as expected. See figure and compare\n');
    fprintf('    to Fig 1 in Newton iter paper by Kressner.\n\n');
    figure; hold on
    plot(phi(linspace(0,2*pi,100)));
    plot(real(evals),imag(evals),'ob'); hold off
    axis([-2.5 0.5 -r r])
    title('Test 9: Example 5.5')
else
    fprintf('  Expected 5 eigenvalues but there are %d\n\n',length(evals));
end


end

function print_results(verbose,B0err,B1err,evalserr,evecserr,k)
  if verbose || norm([B0err,B1err,evalserr,evecserr]) > 1e-6
    fprintf('  k:     %d\n', k);
    fprintf('  B0:    %8.6e\n', B0err);
    fprintf('  B1:    %8.6e\n', B1err);
    fprintf('  evals: %8.6e\n', evalserr);
    fprintf('  evecs: %8.6e\n\n', evecserr);
  end
end