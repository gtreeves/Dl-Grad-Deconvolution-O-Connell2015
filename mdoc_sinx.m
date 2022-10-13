% Set boundary values & step size

clear all
clc
a = 0;
alpha = 0;
b = 3*pi;
beta = 4;

m = 999; % set number of points in f(x)
h = b/(m+1); % set step size

% Create sparse matrix for centered approximation

e = ones(m,1);
A = spdiags([e -2*e e], [-1 0 1], m, m)/h^2;
% Simulate discrete values for f(x) = -sin(x)

x = linspace(a,b,m);
F = -1*sin(x);
F(1) = F(1) - alpha/h^2; % adjust initial value
F(end) = F(end) - beta/h^2; % adjust final value
% Solve system for U

U = A\F';
% Plot U compared to true values

hold on
plot(x,U,'b')

u = sin(x)+x*beta/b;
plot(x,u,'--r')
legend('Approximation', 'True Values','Location','Southeast')

%%
% Set boundary values & step size

clear all
clc
a = 0;
alpha = 0;
b = 4;
beta = 4;

m = 999; % set number of points in f(x)
h = b/(m+1); % set step size
% Create sparse matrix for centered approximation

e = ones(m,1);
A = spdiags([e -2*e e], [-1 0 1], m, m)/h^2;
% Simulate discrete values for f(x) = exp(x)

x = linspace(a,b,m);
F = exp(x);
F(1) = F(1) - alpha/h^2; % adjust initial value
F(end) = F(end) - beta/h^2; % adjust final value
% Solve system for U

U = A\F';

% Plot U compared to true values

figure
hold on
plot(x,U,'b')

u = exp(x)+((5-exp(4))/4)*x-1;
plot(x,u,'--r')
legend('Approximation', 'True Values','Location','Southeast')