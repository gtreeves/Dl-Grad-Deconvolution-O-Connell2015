function [ varargout ] = decon_log_nest_JPattern(Params,varargin)
%DECON_LOG_NEST Nested functions version of Deconvolution
%
%   [ERROR,PENALTY] = decon_log_nest_JPattern(PARAMS) takes a vector of
%   parameters (PARAMS), where size(PARAMS)=[1 20], and calculates the RSS
%   error (ERROR) against those found in GREGSDATA; the PENALTY (= 0 or 1)
%   indicates if the solution is feasible/infeasible, respectively.
%
%   [ERROR,PENALTY,PROTEIN2X,TIME,BETA,NAMES,PROTEIN1X] = ...
%		decon_log_nest_JPattern(PARAMS)
%	solves as above, but creates nested structure PROTEIN2X and PROTEIN1X,
%	in which each molecular species and nuclear cycle is accessible as
%	PROTEIN.SPECIES.NC. The list of species is 'dlNuc', 'dlCyt',
%	'dlCactNuc', 'dlCactCyt', 'cactNuc', and 'cactCyt'. Each species has NC
%	fields 'NC10', 'NC10m', 'NC11', 'NC11m', 'NC12', 'NC12m', 'NC13',
%	'NC13m', and 'NC14'.
%	2X and 1X refer to 2x dl (wt) and 1x dl genotypes. TIME is a similar
%   structure that refers to the accompanying X and T values (spacial
%   coordinates and time vectors) for each species & nuclear cycle. BETA
%   is the value that scales the simulation output to the microscopy data
%   from GREGSDATA. NAMES is a cell variable that contains the names of the
%   parameters in the model.
%

opts = optimset('Display','off');

args = struct('dl', 1, 'cact',1,'plot_type','all');
Beta = 1;
if nargin > 1
    if mod(length(varargin),2)==0
        for i = 1:2:length(varargin)
            args.(varargin{i})=varargin{i+1};
        end
    else
        error('Varargin requires name-value pairs. Some options not specified')
    end
end
options = [];
varargout = cell(1,nargout);
E = zeros(size(Params,1),1); penalty = [E E E];
% Greg's data
[C1,dC1] = gregsData('nuclear');
% misc. initial procedures
ncs = {'NC10','NC10m','NC11','NC11m','NC12','NC12m','NC13','NC13m','NC14'};
if nargout <= 2 && nargout > 0
    names = {'dlNuc','dlCyt','dlCactNuc'};
    times = false;
    
else
    names = {'dlNuc','dlCyt','dlCactNuc','dlCactCyt','cactNuc','cactCyt'};
    ts = {'T10','T10m','T11','T11m','T12','T12m','T13','T13m','T14'};
    xs = {'X10','X10m','X11','X11m','X12','X12m','X13','X13m','X14'};
    nt=1;
    times = true;
end
% options = odeset('RelTol',1e-1,'AbsTol',1e-1);

load('hybrid_embryo.mat','data')
nc = ncs{1};
% Constants
% points in x (# of nuclei)
M = struct('NC10',13,'NC10m',13,'NC11',19,'NC11m',19,'NC12',26,'NC12m',...
    26,'NC13',36,'NC13m',36,'NC14',51);
% Length (L2 - L1) in x
L1 = 0; L2 = 1;
% time spans
tspan = struct('NC10',data.t(1:16),'NC10m',data.t(16:33),'NC11',data.t(1:16)+data.t(33),...
    'NC11m',data.t(17:33)+data.t(33),'NC12',data.t(34:59)+data.t(33),...
    'NC12m',data.t(60:76)+data.t(33),'NC13',data.t(77:126)+data.t(33),...
    'NC13m',data.t(127:148)+data.t(33),'NC14',linspace(35.6,90,184)'+data.t(33));
clear data

%initialize protein
sub_pro = struct('NC10',zeros(13,16),'NC10m',zeros(13,18),'NC11',...
    zeros(19,16),'NC11m',zeros(19,17),'NC12',zeros(26),'NC12m',...
    zeros(26,17),'NC13',zeros(36,50),'NC13m',zeros(36,22),'NC14',...
    zeros(51,184));
protein = struct('dlNuc',sub_pro,'dlCyt',sub_pro,'dlCactNuc',sub_pro,...
    'dlCactCyt',sub_pro,'cactNuc',sub_pro,'cactCyt',sub_pro,...
    'nuclearDorsal',sub_pro);
clear sub_pro

%initialize time
time = struct;

% other nested variables
P=[];An=[];Am=[];Vn=[];Vc=[];


%% MAIN
for rows = 1:size(Params,1)
    
    Params_ = 10.^(Params(rows,:));
    Params_ = num2cell(Params_);
    [lambdaU,lambdaW,lambdaV,sigmaU,sigmaW,sigmaV,muU,muW,muV,...
            gamma,psi,alpha,beta,phi,beta0]=Params_{:};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2x Dorsal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = linspace(0,1,M.NC10)';
    [An,Am,~,Vn,Vc]=nuclearSize(1,'static',M.NC10,'interphase');
    y0 = [zeros(13,1);zeros(13,1);ones(13,1);ones(13,1);
        (1/alpha)*ones(13,1); (1/alpha)*ones(13,1)];
    
    
    for i = 1:9                                                         
        if mod(i,2)==0                                                  
            nc = ncs{i};                                                
            mitosis;                                                    
        else                                                            
            nc = ncs{i};                                                
            interphase;                                                 
        end                                                             
    end                                                                 
    
    % Assign outputs
    if nargout > 2
        
    end
    %combine dlNuc and dlCactNuc
    
    for j=1:length(ncs)
        protein.('nuclearDorsal').(ncs{j})=protein.dlNuc.(ncs{j})+...
            protein.dlCactNuc.(ncs{j});
    end
        protein2x = protein;

    try
        % Simulation data as a column vector
        D1 = [protein.nuclearDorsal.NC11(:);
            protein.nuclearDorsal.NC12(:);
            protein.nuclearDorsal.NC13(:);
            protein.nuclearDorsal.NC14(:)];
        
%         D_free = [protein.dlNuc.NC11(:);
%             protein.dlNuc.NC12(:);
%             protein.dlNuc.NC13(:);
%             protein.dlNuc.NC14(:)];
%         D_bound = [protein.dlCactNuc.NC11(:);
%             protein.dlCactNuc.NC12(:);
%             protein.dlCactNuc.NC13(:);
%             protein.dlCactNuc.NC14(:)];
        
        Beta = mean(D1.*C1)/mean(D1.^2);
        
        E(rows) = sqrt(sum(((C1-Beta*D1)./dC1).^2)/length(D1));
        
        
        dg=protein.nuclearDorsal.NC14(:,end);
        [A,B,C,D]=fitstep(x,dg);
        dgfit = A.*10.^(B)./(10.^(B)+x.^(2*C))+D;
        hmax = min(dgfit)+0.5*(max(dgfit)-min(dgfit));
        sft = @(x) A.*10.^(B)./(10.^(B)+x.^(2*C))+D-hmax;
        zed2=fzero(sft,0.5,opts);
    catch ME
        if strcmp(ME.message,'Matrix dimensions must agree.') || strcmp(ME.message,'Exiting due to infeasibility: 1 lower bound exceeds the corresponding upper bound.')
        E(rows) = 1e6+randn;
        penalty(rows,:)=1e6+randn;
        else
            error(ME.message)
        end
        continue
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1x Dorsal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = linspace(0,1,M.NC10)';
    [An,Am,~,Vn,Vc]=nuclearSize(1,'static',M.NC10,'interphase');
    y0 = [zeros(13,1);zeros(13,1);0.5*ones(13,1);0.5*ones(13,1);
        (1/alpha)*ones(13,1); (1/alpha)*ones(13,1)];
    
    nt = 1;
    for i = 1:9                                                         
        if mod(i,2)==0                                                  
            nc = ncs{i};                                                
            mitosis;                                                    
        else                                                            
            nc = ncs{i};                                                
            interphase;                                                 
        end                                                             
    end                                                                 
    
    % Assign outputs
    if nargout > 2
        
    end
    %combine dlNuc and dlCactNuc
    
    for j=1:length(ncs)
        protein.('nuclearDorsal').(ncs{j})=protein.dlNuc.(ncs{j})+...
            protein.dlCactNuc.(ncs{j});
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
%% Penalty
 
    vtot = protein2x.nuclearDorsal.NC14(1,:); vtot = vtot-min(vtot); vtot = vtot/max(vtot);
    vdl = protein2x.dlNuc.NC14(1,:); vdl = vdl-min(vdl); vdl = vdl/max(vdl);
    
    btot = protein2x.nuclearDorsal.NC14(end,:); btot = btot-min(btot); btot = btot/max(btot);
    bdl = protein2x.dlCactNuc.NC14(end,:); bdl = bdl-min(bdl); bdl = bdl/max(bdl);

    penalty(rows,1) = sqrt((sum((vtot-vdl).^2)+sum((btot-bdl).^2)))-...
        sqrt((sum((vtot-bdl).^2)+sum((btot-vdl).^2)));
    dg=protein.nuclearDorsal.NC14(:,end);
    try
    [A,B,C,D]=fitstep(x,dg);
    dgfit = A.*10.^(B)./(10.^(B)+x.^(2*C))+D;
    hmax = min(dgfit)+0.5*(max(dgfit)-min(dgfit));
    sft = @(x) A.*10.^(B)./(10.^(B)+x.^(2*C))+D-hmax;
    zed1=fzero(sft,0.5,opts);
    z = zed2-zed1;
    catch ME
        z = NaN;
        if ~strcmp(ME.message,'Exiting due to infeasibility: 1 lower bound exceeds the corresponding upper bound.')
            error(ME.message)
        end
    end
    if isnan(z)
        z = 1;
    end
%	penalty(rows,2) = z;
	penalty(rows,2) = penalty(rows,2)+ z+0.03; % adding a small number to force the expansion
	
    if protein2x.dlNuc.NC14(end,end) > protein2x.dlCactNuc.NC14(end,end)
        penalty(rows,3) = 1;
    end
end

%% Output

pnames = {'lambdaU','lambdaW','lambdaV','sigmaU','sigmaW','sigmaV',...
    'muU','muW','muV','gamma','psi','alpha','beta','phi','beta0'};
argsout = {E,penalty,protein2x,time,Beta,pnames,protein};
for i = 1:nargout
    varargout{i} = argsout{i};
end

if nargout == 0
    figure; hold on;
    if strcmp(args.plot_type,'all')
        for i = 1:length(names)
            subplot(2,3,i)
            hold on
            for j = 1:9
            surf(time.(ts{j})',time.(xs{j})',...
                Beta*protein2x.(names{i}).(ncs{j}))
            surf(time.(ts{j})',-time.(xs{j})',...
                Beta*protein2x.(names{i}).(ncs{j}))
            end
            title(names{i})
            shading flat
            view([45 45])
        end
        
    else
        
    for j = 1:9
        surf(time.(ts{j})',time.(xs{j})',...
            Beta*protein2x.(args.plot_type).(ncs{j}))
        surf(time.(ts{j})',-time.(xs{j})',...
            Beta*protein2x.(args.plot_type).(ncs{j}))
    end
    end
    shading flat
    view([45 45])
end

%% nested functions
    function interphase
        %update variables
        x = linspace(L1,L2,M.(nc))';
        e = ones(M.(nc),1);P = spdiags([e -2*e e],[-1 0 1],M.(nc),M.(nc));
        j1 = spdiags(ones(M.(nc),1),0,M.(nc),M.(nc));
        options = odeset('RelTol',1e-1,'AbsTol',1e-2,'JPattern',...
            [P j1 j1 j1 j1 j1 ;
            j1 P j1 j1 j1 j1 ;
            j1 j1 P j1 j1 j1 ;
            j1 j1 j1 P j1 j1 ;
            j1 j1 j1 j1 P j1 ;
            j1 j1 j1 j1 j1 P]);
        
        P(1,2) = 2;P(M.(nc),M.(nc)-1) = 2;
        
        [An,Am,~,Vn,Vc]=nuclearSize(1,'static',M.(nc),'interphase');
        
        % interpolating (nuc = cyt)
        y0_ = reshape(y0,length(y0)/6,6);
        y0 = zeros(M.(nc)*6,1);
        
        for n = 1:2:6
            y0((n-1)*M.(nc)+1:n*M.(nc),1)=interp1(linspace(0,1,...
                size(y0_,1)),y0_(:,n+1),linspace(0,1,M.(nc)));
            y0((n)*M.(nc)+1:(n+1)*M.(nc),1)=interp1(linspace(0,1,...
                size(y0_,1)),y0_(:,n+1),linspace(0,1,M.(nc)));
        end
        
		
        if times && nt <=9
            [t,y] = ode15s(@int_fun,tspan.(nc),y0,options);
            [time.(xs{nt}),time.(ts{nt})] = meshgrid(x,t);
            nt=nt+1;
        else
            [~,y] = ode15s(@int_fun,tspan.(nc),y0,options);
        end
        
        y = y';
		
        
        for i_ = 1:length(names)
            protein.(names{i_}).(nc)=y((i_-1)*M.(nc)+1:(i_)*M.(nc),:);
        end
        
        % update y0 for mitosis (mixing)
        y0 = y(:,end);
        for i_ =[6 5 4 3 2 1]
            if mod(i_,2)==0
                y0((i_-1)*M.(nc)+1:(i_)*M.(nc)) = ...
                    (Vn*y0((i_-2)*M.(nc)+1:(i_-1)*M.(nc)) +...
                    Vc*y0((i_-1)*M.(nc)+1:(i_)*M.(nc)))/(Vn+Vc);
            else
                y0((i_-1)*M.(nc)+1:(i_)*M.(nc))=zeros(M.(nc),1);
            end
		end
        
    end
    function F = int_fun(~,y)
        
    un = y(1:M.(nc));
    uc = y(M.(nc)+1:2*M.(nc));
    wn = y(2*M.(nc)+1:3*M.(nc));
    wc = y(3*M.(nc)+1:4*M.(nc));
    vn = y(4*M.(nc)+1:5*M.(nc));
    vc = y(5*M.(nc)+1:6*M.(nc));
    
    
    
    f1 = (sigmaU*An*uc-muU*An*un-gamma*Vn*(un.*vn)+Vn*beta0*wn)/Vn;
    f2 = (lambdaU*Am*P*uc+Vc*(beta*exp(-0.5*(x/phi).^2)+beta0).*wc-gamma*Vc*(uc.*vc)...
        -sigmaU*An*uc+muU*An*un)/Vc;
    f3 = (sigmaW*An*wc-An*muW*wn+gamma*Vn*(un.*vn)-Vn*beta0*wn)/Vn;
    f4 = (lambdaW*Am*P*wc-Vc*(beta*exp(-0.5*(x/phi).^2)+beta0).*wc+gamma*Vc*(uc.*vc)...
        -sigmaW*An*wc+muW*An*wn)/Vc;
    f5 = (sigmaV*An*vc-muV*An*vn-gamma*psi*Vn*(un.*vn)+psi*Vn*beta0*wn)/Vn;
    f6 = (lambdaV*Am*P*vc+psi*Vc*(beta*exp(-0.5*(x/phi).^2)+beta0).*wc-gamma*psi*Vc*(uc.*vc)+...
        1-alpha*Vc*vc-sigmaV*An*vc+muV*An*vn)/Vc;
    
    F = [f1;f2;f3;f4;f5;f6];
    end

    function mitosis
        % update variables
        [An,Am,~,Vn,Vc]=nuclearSize(1,'static',M.(nc),'mitosis');
        %         options = odeset('Jacobian',@mit_jac);
        if times && nt <=9
            [t,y] = ode15s(@mit_fun,tspan.(nc),y0,options);
            [time.(xs{nt}),time.(ts{nt})] = meshgrid(x,t);
            nt=nt+1;
        else
            [~,y] = ode15s(@mit_fun,tspan.(nc),y0,options);
        end
        y = y';
        
		for i_ = 1:length(names)
			protein.(names{i_}).(nc)=y((i_-1)*M.(nc)+1:(i_)*M.(nc),:);
		end
		
		% update y0 for interphase (interpolation)
		y0 = y(:,end);
        
    end


    function F = mit_fun(~,y)
        
    uc = y(M.(nc)+1:2*M.(nc));
    wc = y(3*M.(nc)+1:4*M.(nc));
    vc = y(5*M.(nc)+1:6*M.(nc));
        
        
        f1 = zeros(M.(nc),1);
        f2 = (lambdaU*Am*P*uc+Vc*(beta*exp(-0.5*(x/phi).^2)+beta0).*wc-gamma*Vc*uc.*vc)/(Vc);
        f3 = zeros(M.(nc),1);
        f4 = (lambdaW*Am*P*wc-Vc*(beta*exp(-0.5*(x/phi).^2)+beta0).*wc+gamma*Vc*uc.*vc)/(Vc);
        f5 = zeros(M.(nc),1);
        f6 = (lambdaV*Am*P*vc+psi*Vc*(beta*exp(-0.5*(x/phi).^2)+beta0).*wc-psi*gamma*Vc*uc.*vc+1-alpha*Vc*vc)/(Vc);
        
        F = [f1;f2;f3;f4;f5;f6];
        
    end
end