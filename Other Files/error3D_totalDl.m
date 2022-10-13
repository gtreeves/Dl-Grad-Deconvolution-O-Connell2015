function error = error3D_totalDl(Model,protein,plotting)

dlNuc = protein{1};
dlCyt = protein{2};
dlCactNuc = protein{3};
dlCactCyt = protein{4};

%% Greg's data
load('hybrid_embryo.mat')
A2 = data.A2;
B2 = data.B2;
t = data.t;
sig = data.Sig2;

% Keeping this as ones, so if we want it later, we still have it in there.
dA = ones(size(A2));
dB = dA;


%% NC 11 (interphase & mitosis)
N11T = t(1:33);
N11A = A2(1:33);
N11B = B2(1:33);
dA11 = dA(1:33);
dB11 = dB(1:33);
sig11 = sig(1:33);

m = -0.07;
M = 19;
x = linspace(0,1,M);
N11D = repmat(N11A,1,M).*exp(-repmat(x,length(N11A),1).^2/2./repmat(sig11,1,M).^2) + ...
    repmat(N11B,1,M) + m*N11A*x;
N11D = fliplr(rot90(N11D,3));

dC11 = exp(-repmat(x,length(N11A),1).^2/2./repmat(sig11,1,M).^2).*repmat(dA11,1,length(x))...
    +m*dA11*x + m*repmat(x,length(dA11),1) + repmat(dB11,1,length(x));

integral11D = trapz(x,N11D);


%% NC 12 (interphase & mitosis)
N12T = t(34:76);
N12A = A2(34:76);
N12B = B2(34:76);
dA12 = dA(34:76);
dB12 = dB(34:76);
sig12 = sig(34:76);
% N12T = t(33:76);
% N12A = A2(33:76);
% N12B = B2(33:76);
% dA12 = dA(33:76);
% dB12 = dB(33:76);
% sig12 = sig(33:76);

M = 26;
x = linspace(0,1,M);
N12D = repmat(N12A,1,M).*exp(-repmat(x,length(N12A),1).^2/2./repmat(sig12...
    ,1,M).^2) + repmat(N12B,1,M) + m*N12A*x;
N12D = fliplr(rot90(N12D,3));
dC12 = exp(-repmat(x,length(N12A),1).^2/2./repmat(sig12,1,M).^2).*...
    repmat(dA12,1,length(x))+m*dA12*x + m*repmat(x,length(dA12),1)...
    + repmat(dB12,1,length(x));

integral12D = trapz(x,N12D);


%% NC 13 (interphase & mitosis)
N13T = t(77:148);
N13A = A2(77:148);
N13B = B2(77:148);
dA13 = dA(77:148);
dB13 = dB(77:148);
sig13 = sig(77:148);
% N13T = t(76:148);
% N13A = A2(76:148);
% N13B = B2(76:148);
% dA13 = dA(76:148);
% dB13 = dB(76:148);
% sig13 = sig(76:148);

M = 36;
x = linspace(0,1,M);
N13D = repmat(N13A,1,M).*exp(-repmat(x,length(N13A),1).^2/2./repmat(sig13,1,M).^2) + ...
    repmat(N13B,1,M) + m*N13A*x;
N13D = fliplr(rot90(N13D,3));
dC13 = exp(-repmat(x,length(N13A),1).^2/2./repmat(sig13,1,M).^2).*repmat(dA13,1,length(x))+m*dA13*x + m*repmat(x,length(dA13),1) + repmat(dB13,1,length(x));

integral13D = trapz(x,N13D);


%% NC 14 (interphase & mitosis)
N14T = t(149:332);
N14A = A2(149:332);
N14B = B2(149:332);
dA14 = dA(149:332);
dB14 = dB(149:332);
sig14 = sig(149:332);

% N14T = t(148:end);
% N14A = A2(148:end);
% N14B = B2(148:end);
% dA14 = dA(148:end);
% dB14 = dB(148:end);
% sig14 = sig(148:end);

M = 51;
x = linspace(0,1,M);
N14D = repmat(N14A,1,M).*exp(-repmat(x,length(N14A),1).^2/2./repmat(sig14,1,M).^2) + ...
    repmat(N14B,1,M) + m*N14A*x;
N14D = fliplr(rot90(N14D,3));
dC14 = exp(-repmat(x,length(N14A),1).^2/2./repmat(sig14,1,M).^2).*...
    repmat(dA14,1,length(x))+m*dA14*x + m*repmat(x,length(dA14),1)...
    + repmat(dB14,1,length(x));

integral14D = trapz(x,N14D);

totalDorsal = [integral11D integral12D integral13D integral14D];

% Fluorescence data as a column vector
C = [N11D(:);N12D(:);N13D(:);N14D(:)];
dC = [dC11(:);dC12(:);dC13(:);dC14(:)];
  


% Combine all dorsal components
Vn10 = 3.43; Vc10 = 44.2; V10 = Vn10+Vc10;
Vn11 = 3.43; Vc11 = 19.0; V11 = Vn11+Vc11;
Vn12 = 2.46; Vc12 = 9.22; V12 = Vn12+Vc12;
Vn13 = 1.37; Vc13 = 4.73; V13 = Vn13+Vc13;
Vn14 = 1;    Vc14 = 2.10; V14 = Vn14+Vc14;

N11S = [(Vc11*dlCactCyt.NC11+Vn11*dlCactNuc.NC11...
    +Vc11*dlCyt.NC11+Vn11*dlNuc.NC11)/V11 dlCyt.NC11m+dlCactCyt.NC11m];
N12S = [(Vc12*dlCactCyt.NC12+Vn12*dlCactNuc.NC12...
    +Vc12*dlCyt.NC12+Vn12*dlNuc.NC12)/V12 dlCyt.NC12m+dlCactCyt.NC12m];
N13S = [(Vc13*dlCactCyt.NC13+Vn13*dlCactNuc.NC13...
    +Vc13*dlCyt.NC13+Vn13*dlNuc.NC13)/V13 dlCyt.NC13m+dlCactCyt.NC13m];
N14S = (Vc14*dlCactCyt.NC14+Vn14*dlCactNuc.NC14...
    +Vc14*dlCyt.NC14+Vn14*dlNuc.NC14)/V14;

integral11S = trapz(linspace(0,1,size(N11S,1)),N11S);
integral12S = trapz(linspace(0,1,size(N12S,1)),N12S);
integral13S = trapz(linspace(0,1,size(N13S,1)),N13S);
integral14S = trapz(linspace(0,1,size(N14S,1)),N14S);

totalDorsalSim = [integral11S integral12S integral13S integral14S];

% Simulation data as a column vector
D = [N11S(:);N12S(:);N13S(:);N14S(:)]; 
    
if (size(D) == size(C))
        Beta = (mean(D.*C)-mean(D)*mean(C))/(mean(D.^2)-mean(D)^2);
        %error = sum((C - Beta*D).^2./dC);
        error = sum(((C - Beta*D)./dC).^2);
    else
        error = 1e12;
    end
        

       if plotting == 1
           figure
           hold on
           plot(C)
           plot(Beta*D)
       end
        if plotting == 2
            figure
            subplot(1,2,1)
            hold on 
            surf(N14T,linspace(0,1,51),N14D)
            surf(N13T,linspace(0,1,36),N13D)
            surf(N12T,linspace(0,1,26),N12D)
            surf(N11T,linspace(0,1,19),N11D)
            shading flat
            
            subplot(1,2,2)
            hold on 
            surf(N14T,linspace(0,1,51),Beta*N14S)
            surf(N13T,linspace(0,1,36),Beta*N13S)
            surf(N12T,linspace(0,1,26),Beta*N12S)
            surf(N11T,linspace(0,1,19),Beta*N11S)
            shading flat
        end

end
