function [ varargout ] = genef_log_nest_avg( Params,dl_grad )
%GENEF_LOG_NEST_AVG Nested version of the gene expression function for the
%deconvolution model, using the average of the last 4 dl gradient time
%points
%   [ERR,PENALTY,PROTEIN] = genef_log_nest(parameters,dorsal_gradient)
if size(Params,1)==15
    Params = Params';
end
if size(Params,2) < 15
    error('GENEF_LOG_NEST requires 15 parameters')
end
% dl_grad = struct;
while iscell(dl_grad)
    dl_grad = dl_grad{1};
end
if ~isstruct(dl_grad)
    error('GENEF_LOG_NEST requires a struct representing the dl gradient.')
end
% misc. initial procedures
ncs = {'NC10','NC10m','NC11','NC11m','NC12','NC12m','NC13','NC13m','NC14'};
names = {'sna','sog','dpp','vnd'};
load('hybrid_embryo.mat','data')
nc = ncs{1};

t = [];  genename = {'sna', 'sog', 'dpp', 'vnd'};
% s0 = [0 0.3 1 0.27]; % from Greg's code
% initialize values for loop
geneavg = zeros(151,4);

for i = 1:4
    matname = [genename{i},'avg.mat'];
    load(matname,'t')
    if i == 2
        geneavg(1:121,i) = t(31:end);
    else
        geneavg(:,i) = t;
    end
    
end
geneavg(1:25,1)=round(geneavg(1:25,1));
geneavg(133:151,3)=round(geneavg(133:151,3));
geneavg = interp1(linspace(0,1,151),geneavg,linspace(0,1,51));

% Constants
% points in x (# of nuclei)
M = struct('NC10',13,'NC10m',13,'NC11',19,'NC11m',19,'NC12',26,'NC12m',...
    26,'NC13',36,'NC13m',36,'NC14',51);
% Length (L2 - L1) in x

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
protein = struct('sna',sub_pro,'sog',sub_pro,'dpp',sub_pro,'vnd',sub_pro,'NC14',[]);
clear sub_pro

% other nested variables
A=[];G=[];D=[];P=[];Vc=[];
E = zeros(size(Params,1),1);
penalty = zeros(size(Params,1),1);
%% MAIN
for rows = 1:size(Params,1)
    
    %     Params_ = 10.^(Params(rows,:));
    %     Params_ = num2cell(Params_);
    Params_ = num2cell(10.^Params(rows,:));
    [thetaDlSna,thetaDlSog,thetaSnaSog,thetaDlDpp,...
        thetaSnaVnd,thetaDlVnd,tauSna,tauSog,tauDpp,tauVnd,...
        nSna,nSog,nDpp,nVnd,noise] = Params_{:};
    noiseInt = 0.05;
    
    sim = zeros(204,10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for gs = 1:10
        [~,~,~,~,Vc]=nuclearSize(1,'static',M.NC10,'interphase');
        y0 = repmat(zeros(M.NC10,1),4,1);
        for i = 1:9                                                         %%
            nc = ncs{i};                                                    %%
            if mod(i,2)==0                                                  %%
                mitosis;                                                    %%
            else                                                            %%
                interphase;                                                 %%
            end                                                             %%
        end                                                                 %%
        sim(:,gs) = [protein.sna.NC14(:,end); protein.sog.NC14(:,end);...
            protein.dpp.NC14(:,end); protein.vnd.NC14(:,end)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sim = mean(sim,2); sim = reshape(sim,51,4);
    % Assign outputs
    if nargout > 2
        
    end
    
    infeasible = false;
    %% Error
    try
        % interpolating down to the size of the simulation data
        
        
        sim(isnan(sim))= 0;
        E(rows) = sum((geneavg(:)-sim(:)).^2);
        protein.NC14 = sim;
    catch
        E(rows) = 1e6+randn;
        infeasible = true;
    end
    % Penalty
    
    if infeasible
        penalty(rows) = 1;
    else
        penalty(rows) = 0;
    end
end
argsout = {E,penalty,protein};
varargout = cell(1,nargout);
for i = 1:nargout
    varargout{i} = argsout{i};
end

%% nested functions
    function interphase
        %update variables
        k = tspan.(nc)(2)-tspan.(nc)(1);
        [~,~,~,~,Vc]=nuclearSize(1,'static',M.(nc),'interphase');
        % interpolating (nuc = cyt)
        y0_ = reshape(y0,length(y0)/4,4);
        y0 = zeros(M.(nc)*4,1);
        
        for n = 1:4
            y0((n-1)*M.(nc)+1:n*M.(nc),1)=interp1(linspace(0,1,...
                size(y0_,1)),y0_(:,n),linspace(0,1,M.(nc)));
        end
        %         options = odeset('Jacobian',@int_jac);
        y = [y0 zeros(length(y0),size(dl_grad.(nc),2)-1)];
        for j = 1:size(dl_grad.(nc),2)-1
            interpDN = repmat(dl_grad.(nc)(:,j),1,4);
            
            noiseDN = interpDN + sqrt(interpDN)*noise.*randn(size(interpDN))+...
                noiseInt.*randn(size(interpDN));
            
            noiseDN(noiseDN<=0)=0.5*min(min(abs(noiseDN)));
            nD=100;
            % Snail Switch - activated by dorsal (noiseDN)
            A = noiseDN(:,1).^nSna./(thetaDlSna^nSna + noiseDN(:,1).^nSna);
            % Sog Switch - activated by dorsal but repressed by
            G = (noiseDN(:,2).^nSog./(thetaDlSog^nSog + noiseDN(:,2).^nSog))...
                .*(thetaSnaSog.^nD./(thetaSnaSog^nD + y(1:M.(nc),j).^nD));
            % Dpp Switch - repressed by dorsal
            P = (thetaDlDpp.^nDpp./(thetaDlDpp^nDpp + noiseDN(:,3).^nDpp));
            % Vnd Switch - activated by dorsal but repressed by sna
            D = (noiseDN(:,4).^nVnd./(thetaDlVnd^nVnd + noiseDN(:,4).^nVnd)).*...
                (thetaSnaVnd.^nD./(thetaSnaVnd^nD + y(1:M.(nc),j).^nD));
            
            y(:,j+1) = [(1/(1+k/tauSna)).*((k/tauSna).*A+y(1:M.(nc),j));
                (1/(1+k/tauSog)).*((k/tauSog).*G+y(M.(nc)+1:2*M.(nc),j));
                (1/(1+k/tauDpp)).*((k/tauDpp).*P+y(2*M.(nc)+1:3*M.(nc),j));
                (1/(1+k/tauVnd)).*((k/tauVnd).*D+y(3*M.(nc)+1:4*M.(nc),j))];
        end
        
        for i_ = 1:length(names)
            protein.(names{i_}).(nc)=y((i_-1)*M.(nc)+1:(i_)*M.(nc),:);
        end
        % update y0 for mitosis
        y0 = y(:,end);
    end


    function mitosis
        % update variables
        [~,~,~,~,Vc]=nuclearSize(1,'static',M.(nc),'mitosis');
        %         options = odeset('Jacobian',@mit_jac);
        
        [~,y] = ode15s(@mit_fun,tspan.(nc),y0);
        
        y = y';
        
        for i_ = 1:length(names)
            protein.(names{i_}).(nc)=y((i_-1)*M.(nc)+1:(i_)*M.(nc),:);
        end
    end

    function y = mit_fun(~,y)
        y = [(-1/tauSna*y(1:M.(nc)))/Vc;
            (-1/tauSog.*y(M.(nc)+1:2*M.(nc)))/Vc;
            (-1/tauDpp.*y(2*M.(nc)+1:3*M.(nc)))/Vc;
            (-1/tauVnd.*y(3*M.(nc)+1:4*M.(nc)))/Vc];
    end
end

