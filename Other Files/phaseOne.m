function [DC,CC,XC,t] = phaseOne(p,Params)

% 
% Unpacking p & Params
% 
[ Tf, Ti]= p{:};

[~,~,~,~,~,~,~,~,~,tauDL,tauC,gamma,xi,kappa,~,~,~]=Params{:};

%
% Create the time span
%
tspan = [Ti Tf];

%
% Create the mesh in x
%

% Initial condition
y0 = [0;0;0];

[t, y] = ode15s(@mitosis,tspan,y0,[],tauDL,tauC,gamma,xi,kappa);

y = y';
DC = y(1,:);
CC = y(2,:);
XC = y(3,:);
 
end


function y = mitosis(t,y,tauDL,tauC,gamma,xi,kappa)
M_cact = 1;
M_dl = dorsalProduction(t);

dc = y(1);
cc = y(2);
xc = y(3);

    f1 = (M_dl-dc-gamma*dc.*cc+kappa*xc)/tauDL;
    f2 = (M_cact-cc-xi*(gamma*dc.*cc-kappa*xc))/tauC;
    f3 = (gamma*dc.*cc-kappa*xc)/tauDL;

    y = [f1;f2;f3];
end