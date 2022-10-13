function C = dcsimple(m,p,Toll)

% gives the concentration profiles for the
% dimensionless version of Kanodia's model for the initial iterate for 
% steadyStateDiffusion3.m given a set of parameters p.

% 
% m = no. of grid points
% p = [dlGamma cactGamma dlCactGamma K1 K2 kappa sigma xi zeta];

%% Unpacking p
 
K1 = p(4);
kappa = p(6); 
sigma = p(7);
xi = p(8); 
zeta = p(9);

%% Inital guess (required by nonlinearity)
dlNuc = 0.5*ones(m,1);
dlCyt = 0.5*ones(m,1);
dlCact = 0.5*ones(m,1);
cactCyt = 0.5*ones(m,1);
% C = [dlNuc;dlCyt;dlCact;cactCyt];
%% Newton's Method

nSteps = 1;  F = 1; dC = 1; tolerance = 1e-10;  nStepsMax = 20;
while norm([F;dC]) > tolerance &&  nSteps <= nStepsMax
	
	%
	% Putting together F
	%
	f1 = dlCyt - K1*dlNuc; 
    f2 = Toll'.*dlCact - kappa*dlCyt.*cactCyt - sigma*(dlCyt - K1*dlNuc);
    f3 = -Toll'.*dlCact + kappa*dlCyt.*cactCyt; 
    f4 = Toll'.*dlCact - kappa*dlCyt.*cactCyt - xi*cactCyt + zeta;
    F = [f1;f2;f3;f4];
    
	%
    % The Jacobian
	%
	I = speye(m); Z = sparse(m,m);
    j1 = [-K1*I  I  Z  Z];
	
	j21 = sigma*K1*I;
	j22 = -kappa*spdiags(cactCyt,0,m,m)-sigma*I;
	j23 = spdiags(Toll',0,m,m);
	j24 = -kappa*spdiags(dlCyt,0,m,m);
    j2 = [j21 j22 j23 j24];
	
	j32 = kappa*spdiags(cactCyt,0,m,m);
	j33 = -spdiags(Toll',0,m,m);
	j34 = kappa*spdiags(dlCyt,0,m,m);
    j3 = [Z j32 j33 j34];
	
	j42 = -kappa*spdiags(cactCyt,0,m,m);
	j43 = spdiags(Toll',0,m,m);
	j44 = -kappa*spdiags(dlCyt,0,m,m)-xi*I;
    j4 = [Z j42 j43 j44];
    J = [j1;j2;j3;j4];
    
	%
    % Step
	%
    dC = -J\F;
	
	dlNuc = dlNuc + dC(1:(m),1);
	dlCyt = dlCyt + dC(m+1:2*(m),1);
	dlCact = dlCact + dC(2*(m)+1:3*(m),1);
	cactCyt = cactCyt + dC(3*(m)+1:4*(m),1);
	
	nSteps = nSteps + 1;


end

C = [dlNuc;dlCyt;dlCact;cactCyt];
end