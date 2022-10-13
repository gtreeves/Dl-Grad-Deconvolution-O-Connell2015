function [A,B,C,D] = fitstep_4ps(s,t)

% FITSTEP [A,B,C,D] = fitstep(x,y)
% t is intensity, s is position in space

[tmax,~] = max(t); tmin = min(t);

%
% The initial guesses and upper and lower bounds of each parameter
%
A = tmax; 
B = -2; 
C = 6; 
D = tmin;
if min(A) < 2e-6 
    A = 0;
    B = 0;
    return
end

AL = 0.1*A; AU = 10*A;
BL = -100; BU = 100;
CL = 1; CU = 100;
DL = 0.1*D; DU = 10*D;
%
% Defining options and the equation to fit to
%
f = fittype('A.*10.^(B)./(10.^(B)+x.^(2*C))+D');
opts = fitoptions('Method','NonlinearLeastSquares',...
	'Startpoint',[A B C D],...
	'Lower',[AL BL CL DL],...
	'Upper',[AU BU CU DU]);

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
C = coeffvals(3);
D = coeffvals(4);

% t_ = A.*10.^(B)./(10.^(B)+s.^(2*C))+D;
% figure
% plot(s,t,'-b',s,t_,'-.r')


