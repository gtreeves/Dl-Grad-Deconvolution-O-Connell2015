%% fullKanodia.m -- This script solves the no-diffusion condition in
% time and space.

% n = dlNuc
% u = dlCyt
% w = dlCact
% v = cactCyt

%
% Mesh in x and t
%
M = 20; % points in x
Tau = 80; % points in t
L = 1; % Length in x
T = 8; % Time span

x = linspace(0,L,M)';
h = x(2)-x(1);

t = linspace(0,T,Tau)';
dT = t(2)-t(1);

%
% Adjustable Parameters 
%
sigma = 0.1; mu = 0.1; lambda = 0.01; gamma = 1; alpha = 30; psi = 0.5; 
beta = 0.0001; phi = 0.0001; xi = 8;

%
% Toll signaling 
%
fx = beta./(phi+x.^xi);

%
%Initial conditions/matrix initialization
%
N=zeros(M,Tau); U = N; W = ones(M,Tau); V = W;

%
% Necessary matrices
%
I = speye(M); Z = spdiags(zeros(M,1),0,M,M);
e = ones(M,1); P = spdiags([e -2*e e],[-1 0 1],M,M); P(1,2) = 2; P(M,M-1) = 2;
P = P/h^2;

for j=1:(Tau-1)
    
    %
    % Updated initial guess for Newton's method
    %
    n = N(:,j+1); n0 = N(:,j);
    u = U(:,j+1); u0 = U(:,j);
    w = W(:,j+1); w0 = W(:,j);
    v = V(:,j+1); v0 = V(:,j);
    
    %
    % Initialize values
    %
    F = 1; del = 1; tolerance = 1e-6; nSteps = 1; nStepsMax = 200;
    while norm([F;del]) > tolerance && nSteps < nStepsMax
%
% Function
%
    %
    % Operations to speed up the code
    %
    tollF = fx.*w;
    tollF0 = fx.*w0;

    %
    % Construct F
    %
    f1 = 0.5*(sigma*u-mu*n+sigma*u0-mu*n0)-(n-n0)/dT;
    f2 = 0.5*(lambda*P*u+tollF-gamma*u.*v-sigma*u+mu*n+...
        lambda*P*u0+tollF0-gamma*u0.*v0-sigma*u0+mu*n0)-(u-u0)/dT;
    f3 = 0.5*(lambda*P*w-tollF+gamma*u.*v+lambda*P*w0-tollF0+gamma*u0.*v0)-(w-w0)/dT;
    f4 = 0.5*(gamma*P*v+psi*tollF-gamma*psi*u.*v-alpha*v...
        +gamma*P*v0+psi*tollF0-gamma*psi*u0.*v0-alpha*v0)+1-(v-v0)/dT;
    F = [f1;f2;f3;f4];

% 
% Jacobian
%
    %
    % Repetative operations that slow down the code
    %
    tollJ = spdiags(fx,0,M,M);
    sigmaJ = sigma*I;
    uJ = spdiags(u,0,M,M);
    vJ = spdiags(v,0,M,M);
    
    %
    % Construct the Jacobian
    %
    j1 = [ -(mu/2+1/dT)*I 0.5*sigmaJ Z Z];
    j2 = [ 0.5*sigmaJ  0.5*(lambda*P-gamma*vJ-sigmaJ)-I/dT ...
        0.5*tollJ -gamma*0.5*uJ];
    j3 = [ Z gamma*0.5*vJ 0.5*(gamma*P-tollJ)-I/dT ...
          0.5*gamma*uJ];
    j4 = [ Z -0.5*gamma*psi*vJ 0.5*psi*tollJ ...
         0.5*(gamma*P-gamma*psi*uJ-alpha*I)-I/dT];
    J = [j1;j2;j3;j4];
    
    del = -J\F;
    deltaN = del(1:(M),1); newN = n + deltaN;
    deltaU = del(M+1:2*(M),1); newU = u + deltaU;
    deltaW = del(2*(M)+1:3*(M),1); newW = w + deltaW;
    deltaV = del(3*(M)+1:4*(M),1); newV = v + deltaV;
       
        %
        % To help prevent negative values...
        %
        if any([newN;newU;newW;newV] < 0) 
            deltaN = deltaN/2;
            newN = n + deltaN;
          
            deltaU = deltaU/2;
            newU = u + deltaU;
        
            deltaW = deltaW/2;
            newW = w + deltaW;
        
            deltaV = deltaV/2;
            newV = v + deltaV;
        end
        
    %
    % Update values
    %
    N(:,j+1) = newN;
    U(:,j+1) = newU;
    W(:,j+1) = newW;
    V(:,j+1) = newV;
   nSteps = nSteps + 1;
    end

end
subplot(2,2,1)
surf(N)
hold on
ylabel('space')
xlabel('time')
title('dlNuc')
shading flat

subplot(2,2,2)
surf(U)
hold on
ylabel('space')
xlabel('time')
title('dlCyt')
shading flat

subplot(2,2,3)
surf(W)
hold on
ylabel('space')
xlabel('time')
title('dlCact')
shading flat

subplot(2,2,4)
surf(V)
hold on
ylabel('space')
xlabel('time')
title('cactCyt')
shading flat

