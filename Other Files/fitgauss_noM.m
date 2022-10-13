function [A,B,M,sig] = fitgauss_noM(s,t)

% t is intensity, s is position in space

[tmax,~] = max(t); tmin = min(t);

%
% The initial guesses and upper and lower bounds of each parameter
%
A = tmax - tmin; 
B = tmin; 
M = s(find(t==max(t),1,'first'));
sig = range(s)/6; sigL = sig/10; sigU = 10*sig;

if min(A) < 2e-6 
    A = 0;
    B = 0;
    M = 0;
    sig = 0;
    return
end

AL = 0.1*A; AU = 10*A;
BL = 0; BU = A + tmin;
ML = min(s); MU = max(s);
%
% Defining options and the equation to fit to
%
% f = fittype('A*exp(-(x)^2/2/sigma^2)+B+M*abs(x)');
f = fittype('A*exp(-(x-M)^2/2/sigma^2)+B');
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
[cfun,~] = fit(s,t,f,opts);

%
% Parameter values
%
coeffvals = coeffvalues(cfun);
A = coeffvals(1);
B = coeffvals(2);
M = coeffvals(3);
sig = coeffvals(4);
