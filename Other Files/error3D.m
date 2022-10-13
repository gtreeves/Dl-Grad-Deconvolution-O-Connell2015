function error = totalError(protein,plotting,w)

dlNuc = protein{1};
dlCyt = protein{2};
dlCactNuc = protein{3};
dlCactCyt = protein{4};

%% Nuclear dorsal

% Greg's data
[C1,dC1] = gregsData('nuclear');

% Combine simulation data
N11S = dlNuc.NC11 + dlCactNuc.NC11;
N12S = dlNuc.NC12 + dlCactNuc.NC12;
N13S = dlNuc.NC13 + dlCactNuc.NC13;
N14S = dlNuc.NC14 + dlCactNuc.NC14;

% weighting matrix
w11 = [ones(9,16); w*ones(10,16)];
w12 = [ones(13,26); w*ones(13,26)];
w13 = [ones(18,50); w*ones(18,50)];
w14 = [ones(25,184); w*ones(26,184)];
w = [w11(:);w12(:);w13(:);w14(:)];

% Simulation data as a column vector
D1 = [N11S(:);N12S(:);N13S(:);N14S(:)];

%% Total Dorsal
% Greg's data
[C2,dC2] = gregsData('nuclear');

% Combine all dorsal components
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

% integral11S = trapz(linspace(0,1,size(N11S,1)),N11S);
% integral12S = trapz(linspace(0,1,size(N12S,1)),N12S);
% integral13S = trapz(linspace(0,1,size(N13S,1)),N13S);
% integral14S = trapz(linspace(0,1,size(N14S,1)),N14S);
% 
% totalDorsalSim = [integral11S integral12S integral13S integral14S];

% Simulation data as a column vector
D2 = [N11S(:);N12S(:);N13S(:);N14S(:)]; 

%% Combine nuclear and total Dorsal column vectors

D = [D1;D2];
C = [C1;C2];
dC = [dC1;dC2];

l = length(w);
W = ones(size(dC)); 
W(1:l) = w;

%% Error calculation
if (size(D) == size(C))
    Beta = abs((mean(D.*C)-mean(D)*mean(C))/(mean(D.^2)-mean(D)^2));
    error = sum(((C - Beta*D).*W./dC).^2);
else
    error = 1e12;
end

%% Optional plotting
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
