function [C] = steadyStateNoDiffusionf(m,p,Toll)

% steadyStateNoDiffusion gives the concentration profiles for the
% dimensionless version of Kanodia's model for the initial iterate for 
% steadyStateDiffusion3.m given a set of parameters p.

% 
% m = no. of grid points
% p = [dlGamma cactGamma dlCactGamma K1 K2 kappa sigma xi zeta];

%% Unpacking p
 
K1 = p(4);
kappa = p(6); 
xi = p(8); 
zeta = p(9);

%% Inital guess (required by nonlinearity)
dlNuc = 0.5*ones(m,1);
dlCyt = 0.5*ones(m,1);
dlCact = 0.5*ones(m,1);
cactCyt = 0.5*ones(m,1);
C = [dlNuc;dlCyt;dlCact;cactCyt];
%% Newton's Method

nSteps = 1;  F = 1; dC = 1; tolerance = 1e-10;  nStepsMax = 20;
while norm([F;dC]) > tolerance &&  nSteps <= nStepsMax
% Values
f1 = dlCyt - K1*dlNuc;
f2 = Toll'.*dlCact - kappa*dlCyt.*cactCyt;
f3 = xi*cactCyt - zeta;
f4 = dlCyt + dlNuc + dlCact - 1;
F = [f1;f2;f3;f4];

% Jacobian
I = speye(m);
Z = spdiags(zeros(m,1),0,m,m);

j1 = [ -spdiags(K1*ones(m,1),0,m,m) I Z Z];
j2 = [ Z spdiags(-kappa*cactCyt,0,m,m) spdiags(Toll',0,m,m) spdiags(-kappa*dlCyt,0,m,m)];
j3 = [ Z Z Z spdiags(xi*ones(m,1),0,m,m)];
j4 = [ I I I Z ];
J = [j1;j2;j3;j4];

dC = -J\F;

dlNuc = dlNuc + dC(1:(m),1);
dlCyt = dlCyt + dC(m+1:2*(m),1);
dlCact = dlCact + dC(2*(m)+1:3*(m),1);
cactCyt = cactCyt + dC(3*(m)+1:4*(m),1);
C = [dlNuc;dlCyt;dlCact;cactCyt];

nSteps = nSteps + 1;


end

end