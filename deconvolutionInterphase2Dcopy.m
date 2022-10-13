function [UN,UC,WN,WC,VN,VC,t] = deconvolutionInterphase2Dcopy

tspan=0:90;
nuclei='static';
% Params={ 0.0210 0.0590 0.0552 0.7792 0.1087 0.1592 0.0666 0.1818 4.4358...
% 0.0570 68.5070 0.0249 1.0356 0.1153 0.0108};
Params={  1 10 10 0.7792 0.1087 0.1592 0.0666 0.1818 4.4358 0.0570...
    68.5070 0.0249 1.0356 0.1153 0.0108};


M=25;
N=floor(M*1.39);
% M = 50; N = 25;
L1=0;
L2=1;
An = 1; Am = 3.3716; Vn = 1; Vc = 14.6746; dVndt = 0;



% Create the mesh in x
%
x = linspace(L1,L2,N)';


h=1;
g=1;
e1 = ones(1,M); 
e2 = ones(1,N); 

e1(end) = 0; e1(end-1) = 2; e1 = repmat(e1,N,1); e1 = e1(:);
e2(end) = 0; e2(end-1) = 2; e2 = repmat(e2',M,1); e2 = e2(:);
e3 = ones(M*N,1);
e4 = flipud(e2);
e5 = flipud(e1);
P = spdiags([e1/g^2 e2/h^2 -2*e3/g^2-2*e3/h^2 e4/h^2 e5/g^2],[-N -1 0 1 N],M*N,M*N);
clear e1 e2 e3 e4 e5


%
% Initial Condition
%
UN = zeros(M*N,1);
UC = UN;
WN = ones(M*N,1);

y0 = [UN(:,1);UC(:,1);WN(:,1);WN(:,1);WN(:,1);WN(:,1)];



tic %'RelTol',1e-4,'AbsTol',1e-6,
options = odeset('RelTol',1e-4,'AbsTol',1e-6,'JPattern',jacpat(M,N));
% try
    [t,y] = ode15s(@interphase,tspan,y0,options,Params,tspan(1),tspan(end),...
        nuclei,M,P,repmat(x,M,1),An,Am,Vn,Vc,dVndt,N);
% catch
%     UN = NaN(size(UN)); UC = UN; WN = UN; WC = UN; VN = UN; VC = UN;
%     t = 0;
%     warning('integration took too long')
%     return
% end
y = y';

UN = y(1:M*N,:);
UC = y(M*N+1:2*M*N,:);
WN = y(2*M*N+1:3*M*N,:);
WC = y(3*M*N+1:4*M*N,:);
VN = y(4*M*N+1:5*M*N,:);
VC = y(5*M*N+1:6*M*N,:);

surf_vect = UN(:,end);
surf_rect = reshape(surf_vect,N,M);
figure
surf(linspace(0,1,M),linspace(0,1,N),surf_rect)
shading flat
hold on
v = [1.05 1.1 1.15]*min(surf_rect(1,:));
contour3(linspace(0,1,M),linspace(0,1,N),surf_rect,[v,v],'k')%,'Color','black','LineWidth',3)
set(gca,'View',[0 90])
end

function F = interphase(t,y,Params,~,~,nuclei,M,P,x,An,Am,Vn,Vc,dVndt,N)
% timeElapsed = toc;
% if timeElapsed > 20
%     
% else
    if length(Params)>14
        [lambdaU,lambdaW,lambdaV,sigmaU,sigmaW,sigmaV,muU,muW,muV,...
            gamma,psi,alpha,beta,phi,beta0]=Params{:};
    else
        [lambdaU,lambdaW,lambdaV,sigmaU,sigmaW,sigmaV,muU,muW,muV,...
            gamma,psi,alpha,beta,phi]=Params{:};
        beta0=0;
    end
    if strcmp(nuclei,'dynamic')
        [An,Am,~,Vn,Vc,dVndt] = nuclearSize(t,nuclei,M,'interphase2D');
    end
    
    dVcdt = -dVndt;
    
    un = y(1:M*N,:);
    uc = y(M*N+1:2*M*N,:);
    wn = y(2*M*N+1:3*M*N,:);
    wc = y(3*M*N+1:4*M*N,:);
    vn = y(4*M*N+1:5*M*N,:);
    vc = y(5*M*N+1:6*M*N,:);
    
    x0 = repmat(linspace(0,1,M),N,1);
    x0= x0(:);
%     x=x-0.5;
%     x0=x0-0.5;
    
    
    tol=getToll(x,x0,phi,2);
    f1 = (sigmaU*An*uc-muU*An*un-gamma*Vn*(un.*vn)+Vn*beta0*wn-un*dVndt)/Vn;
    f2 = (lambdaU*Am*P*uc+Vc*(beta*tol+beta0).*wc...
        -gamma*Vc*(uc.*vc)-sigmaU*An*uc+muU*An*un-uc*dVcdt)/Vc;
    f3 = (sigmaW*An*wc-An*muW*wn+gamma*Vn*(un.*vn)-Vn*beta0*wn-wn*dVndt)/Vn;
    f4 = (lambdaW*Am*P*wc-Vc*(beta*tol+beta0).*wc...
        +gamma*Vc*(uc.*vc)-sigmaW*An*wc+muW*An*wn-wc*dVcdt)/Vc;
    f5 = (sigmaV*An*vc-muV*An*vn-gamma*psi*Vn*(un.*vn)+psi*Vn*beta0*wn-vn*dVndt)/Vn;
    f6 = (lambdaV*Am*P*vc+psi*Vc*(beta*tol+beta0).*wc...
        -gamma*psi*Vc*(uc.*vc)+...
        1-alpha*Vc*vc-sigmaV*An*vc+muV*An*vn-vc*dVcdt)/Vc;
    
    F = [f1;f2;f3;f4;f5;f6];
% end
% if M ==51
% scatter(t,f3(1))
% pause(0.1)
% end
end

function pat = jacpat(M,N)
e1 = ones(1,M); e1 = repmat(e1,N,1); e1 = e1(:);
e2 = ones(1,N); e2 = repmat(e2',M,1); e2 = e2(:);
e3 = ones(M*N,1);
e4 = flipud(e2);
e5 = flipud(e1);
P = spdiags([e1 e2 e3 e4 e5],[-N -1 0 1 N],M*N,M*N);
D = speye(M*N);
Z = 0.*D;

pat = [D D D Z D Z
    D P Z D Z D
    D Z D D D Z
    Z D D P Z D
    D Z D Z D D
    Z D Z D D P];
end

function toll = getToll(x,y,phi,n)
    toll = exp(-0.5*((x)/phi).^n)+2*exp(-0.5*((y)/(phi)).^n);
%     toll = 2*exp(-0.5*(y/phi).^n);
end