clk = clock;
id = strcat('results',num2str(clk(6)),'decon_log_nest_JPattern');

fhandle = 'decon_log_nest_JPattern'; mm = 'min';
lambda = 175; mu = ceil(lambda/7);
pf = .33; varphi = 1; G = 1000; tmax=24*3600;

% [lambdaU,lambdaW,lambdaV,sigmaU,sigmaW,sigmaV,muU,muW,muV,...
%            gamma,psi,alpha,beta,phi,beta0]=Params_{:};
lu = [-5*ones(1,15);5*ones(1,15)];
lu(1,14)=-1;
% lu = [-3*ones(1,20); 2*ones(1,20)]; 
% lu(1,14)=-1; lu(2,14)=-.699; % lower and upper limits

[xb,Statistics,Gm]...
    =isres(fhandle,mm,lu,lambda,G,mu,pf,varphi,tmax);

try
    save(strcat('/share2/mdoconne/deconvolutionResults/',id,'.mat'))
catch
    save(strcat(id,'.mat'))
end
