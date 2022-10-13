%% steadyStateDiffusion3.m -- Nondimensional version of steadyStateDiffusion2.m

% This program simulates the Kanodia model at steady state without neglecting
% diffusion. The parameters have been assigned arbitrary values having
% nothing to do with the actual kinetics of the system. The nonlinearity
% requires the use of Newton's method to reach a solution. 
clear all
clc

%% Grid spacing
m = 50;
h = 1/(m+1);
x = linspace(0,1,m); % m gives 2 extra lines for boundary condition equations
e = ones(1,m);

%% Dimensionless parameters
Toll_i = 0.000001 + x.^4; Toll = 1./Toll_i; % Creates f(x), the dimensionless Toll signal profile
dlGamma = 0.01; K1 = 1.1; K2 = 0.0003; kappa = 1; 
sigma = 1; xi = 5; zeta = 1; cactGamma = 0.01; dlCactGamma = 0.01;
% Vector of parameters for passing to a function
p = [dlGamma; cactGamma; dlCactGamma; K1; K2; kappa; sigma; xi; zeta]; 

%% Initial Iterate
[C] = steadyStateNoDiffusionf_gtr(m,p,Toll);
dlNuc = C(1:(m),1);
dlCyt = C(m+1:2*(m),1);
dlCact = C(2*(m)+1:3*(m),1);
cactCyt = C(3*(m)+1:4*(m),1);

%% Plot initial iterate (steady state, no diffusion solution):
subplot(2,1,1)
hold on
plot(x,[dlNuc dlCyt dlCact cactCyt]);
legend('dlNuc','dlCyt','dlCact','cactCyt','Location','NorthEast')
title('Steady-State, No Diffusion')


%% Run Newton's method for steady state WITH diffusion, using initial iterate as above:
nSteps = 1; G = 1; dC = 1; tolerance = 1e-8; nStepsMax = 50;
e = ones(m,1); P = spdiags([e -2*e e],[-1 0 1],m,m); 
P(1,2) = 2; P(m,m-1) = 2; % Create the sparse tribanded matrix
I = speye(m); Z = sparse(m,m); % Create the Identity matrix and a sparse 
% diagonal matrix of zeros
while norm([G;dC]) > tolerance && nSteps < nStepsMax
    % Matrix G of values 
    g1 = dlCyt - K1*dlNuc; 
    g2 = dlGamma*P*dlCyt/h^2 + Toll'.*dlCact - kappa*dlCyt.*cactCyt - ...
		sigma*(dlCyt - K1*dlNuc);
    g3 = dlCactGamma*P*dlCact/h^2 - Toll'.*dlCact + kappa*dlCyt.*cactCyt; 
    g4 = cactGamma*P*cactCyt/h^2 + Toll'.*dlCact - kappa*dlCyt.*cactCyt - ...
		xi*cactCyt + zeta;
    G = [g1;g2;g3;g4];
    
    % Jacobian J
    j1 = [-K1*I  I  Z  Z];
	
	j21 = sigma*K1*I;
	j22 = dlGamma*P/h^2-kappa*spdiags(cactCyt,0,m,m)-sigma*I;
	j23 = spdiags(Toll',0,m,m);
	j24 = -kappa*spdiags(dlCyt,0,m,m);
    j2 = [j21 j22 j23 j24];
	
	j32 = kappa*spdiags(cactCyt,0,m,m);
	j33 = dlCactGamma*P/h^2-spdiags(Toll',0,m,m);
	j34 = kappa*spdiags(dlCyt,0,m,m);
    j3 = [Z j32 j33 j34];
	
	j42 = -kappa*spdiags(cactCyt,0,m,m);
	j43 = spdiags(Toll',0,m,m);
	j44 = cactGamma*P/h^2-kappa*spdiags(dlCyt,0,m,m)-xi*I;
    j4 = [Z j42 j43 j44];
    J = [j1;j2;j3;j4];
    
    % Step
    dC = -J\G;
    dlNuc = dlNuc + dC(1:(m),1);
    dlCyt = dlCyt + dC(m+1:2*(m),1);
    dlCact = dlCact + dC(2*(m)+1:3*(m),1);
    cactCyt = cactCyt + dC(3*(m)+1:4*(m),1);
    
    nSteps = nSteps + 1;
    
end

subplot(2,1,2)
plot(x,[dlNuc dlCyt dlCact cactCyt]);
legend('dlNuc','dlCyt','dlCact','cactCyt','Location','NorthEast')
title('Steady-State, Diffusion')

