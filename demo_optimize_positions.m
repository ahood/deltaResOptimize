% This demo shows how to use optimize_positions.m. Be aware that poor
% initial position/k0 combinations yield poor results.

%% Choose initial configuration

L = [1,1];
p = [-1,1];

%% Look at where the resonances are

ks = resonances_chebsol(L,p);
figure, plot(real(ks),imag(ks),'*');
axis([-50 50 -5 1])

%% Choose k0 so that the Newton iteration is likely to work

%%
% It looks like there's a resonance pretty near 10 - 1.5i, so 
% using that value for k0 seems reasonable.
axis([9 11 -2.5 -0.5])
k0 = 10-1.5i;

%%
axis([0 20 -5 1])
hold on, plot(real(k0),imag(k0),'ro');

%% See if optimize_positions.m got us closer to k0

[alpha,k] = optimize_positions(k0,L,p);
ks_final = resonances_chebsol(L,alpha);
hold on, plot(real(ks_final),imag(ks_final),'+g');

%%
% Yes, we got closer
axis([9 11 -2.5 -0.5])