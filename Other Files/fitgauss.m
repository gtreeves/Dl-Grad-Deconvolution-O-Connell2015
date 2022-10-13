function [A,B,M,sig] = fitgauss(s,t)

% FITGAUSS [A,B,M,sig] = fitgauss(s,t)
% t is intensity, s is position in space

[tmax,~] = max(t); tmin = min(t);

%
% The initial guesses and upper and lower bounds of each parameter
%
A = tmax - tmin; 
B = tmin; 
M = 0; ML = -1e6; MU = 1e6;
sig = 0.15; sigL = 0.01; sigU = 1;

if min(A) < 2e-6 
    A = 0;
    B = 0;
    M = 0;
    sig = 0;
    return
end

AL = 0.1*A; AU = 10*A;
BL = 0; BU = A + tmin;
%
% Defining options and the equation to fit to
%
f = fittype('A*exp(-(x)^2/2/sigma^2)+B+M*abs(x)');
% f = fittype('A*exp(-(x)^2/2/sigma^2)+B');
opts = fitoptions('Method','NonlinearLeastSquares',...
	'Startpoint',[A B M sig],...
	'Lower',[AL BL ML sigL],...
	'Upper',[AU BU MU sigU]);
% opts = fitoptions('Method','NonlinearLeastSquares',...
% 	'Startpoint',[A B sig],...
% 	'Lower',[AL BL sigL],...
% 	'Upper',[AU BU sigU]);
%
% The actual fit
%
if size(s,1) < size(s,2)
    s=s';
end
if size(t,1) < size(t,2)
    t=t';
end
[cfun,~] = fit(s,t,f,opts);

%
% Parameter values
%
coeffvals = coeffvalues(cfun);
A = coeffvals(1);
B = coeffvals(2);
M = coeffvals(3);
sig = coeffvals(4);

