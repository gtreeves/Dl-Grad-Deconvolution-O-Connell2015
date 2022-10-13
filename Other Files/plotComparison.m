function error = plotComparison(dlNuc,dlCactNuc,time)

%% Greg's data
load('dorsalgradient.mat')
A = data.A;
B = data.B;
t = data.t;
%% NC 11
N11T = t(1:16);
N11A = A(1:16);
N11B = B(1:16);
sig = 0.15;
m = -0.07;
M = 17;
x = linspace(0,1,M);
N11D = repmat(N11A,1,M).*exp(-repmat(x,length(N11A),1).^2/2/sig^2) + ...
    repmat(N11B,1,M) + m*N11A*x;
N11D = fliplr(rot90(N11D,3));

%% NC 12
N12T = t(33:60);
N12A = A(33:60);
N12B = B(33:60);
M = 25;
x = linspace(0,1,M);
N12D = repmat(N12A,1,M).*exp(-repmat(x,length(N12A),1).^2/2/sig^2) + ...
    repmat(N12B,1,M) + m*N12A*x;
N12D = fliplr(rot90(N12D,3));

%% NC 13
N13T = t(76:126);
N13A = A(76:126);
N13B = B(76:126);
M = 37;
x = linspace(0,1,M);
N13D = repmat(N13A,1,M).*exp(-repmat(x,length(N13A),1).^2/2/sig^2) + ...
    repmat(N13B,1,M) + m*N13A*x;
N13D = fliplr(rot90(N13D,3));

%% NC 14
N14T = t(147:327);
N14A = A(147:327);
N14B = B(147:327);
M = 51;
x = linspace(0,1,M);
N14D = repmat(N14A,1,M).*exp(-repmat(x,length(N14A),1).^2/2/sig^2) + ...
    repmat(N14B,1,M) + m*N14A*x;
N14D = fliplr(rot90(N14D,3));


%% Error based on dorsal/ventral midline data & normalized gradients:
    %
    % Midlines:
    %
    totalDlMid1 = dlNuc.NC11(1,:)+dlCactNuc.NC11(1,:);
    totalDlMid2 = dlNuc.NC12(1,:)+dlCactNuc.NC12(1,:);
    totalDlMid3 = dlNuc.NC13(1,:)+dlCactNuc.NC13(1,:);
    totalDlMid4 = dlNuc.NC14(1,:)+dlCactNuc.NC14(1,:);
    
    totalDlBasal1 = dlNuc.NC11(end,:)+dlCactNuc.NC11(end,:);
    totalDlBasal2 = dlNuc.NC12(end,:)+dlCactNuc.NC12(end,:);
    totalDlBasal3 = dlNuc.NC13(end,:)+dlCactNuc.NC13(end,:);
    totalDlBasal4 = dlNuc.NC14(end,:)+dlCactNuc.NC14(end,:);
    
    venusMid1 = N11D(1,:);
    venusMid2 = N12D(1,:);
    venusMid3 = N13D(1,:); 
    venusMid4 = N14D(1,:);
    venusBasal1 = N11D(end,:);
    venusBasal2 = N12D(end,:);
    venusBasal3 = N13D(end,:);
    venusBasal4 = N14D(end,:);

    totalDlMid1 = interp1((time.T11(:,1)-min(time.T11(:,1)))/...
        (max((time.T11(:,1)-min(time.T11(:,1))))),totalDlMid1,...
        (N11T-min(N11T))/(max(N11T-min(N11T))));
    totalDlMid2 = interp1((time.T12(:,1)-min(time.T12(:,1)))/...
        (max((time.T12(:,1)-min(time.T12(:,1))))),totalDlMid2,...
        (N12T-min(N12T))/(max(N12T-min(N12T))));
    totalDlMid3 = interp1((time.T13(:,1)-min(time.T13(:,1)))/...
        (max((time.T13(:,1)-min(time.T13(:,1))))),totalDlMid3,...
        (N13T-min(N13T))/(max(N13T-min(N13T))));
    totalDlMid4 = interp1((time.T14(:,1)-min(time.T14(:,1)))/...
        (max((time.T14(:,1)-min(time.T14(:,1))))),totalDlMid4,...
        (N14T-min(N14T))/(max(N14T-min(N14T))));
    
    totalDlBasal1 = interp1((time.T11(:,1)-min(time.T11(:,1)))/...
        (max((time.T11(:,1)-min(time.T11(:,1))))),totalDlBasal1,...
        (N11T-min(N11T))/(max(N11T-min(N11T))));
    totalDlBasal2 = interp1((time.T12(:,1)-min(time.T12(:,1)))/...
        (max((time.T12(:,1)-min(time.T12(:,1))))),totalDlBasal2,...
        (N12T-min(N12T))/(max(N12T-min(N12T))));
    totalDlBasal3 = interp1((time.T13(:,1)-min(time.T13(:,1)))/...
        (max((time.T13(:,1)-min(time.T13(:,1))))),totalDlBasal3,...
        (N13T-min(N13T))/(max(N13T-min(N13T))));
    totalDlBasal4 = interp1((time.T14(:,1)-min(time.T14(:,1)))/...
        (max((time.T14(:,1)-min(time.T14(:,1))))),totalDlBasal4,...
        (N14T-min(N14T))/(max(N14T-min(N14T))));
   
     totalDlMid = [totalDlMid1; totalDlMid2; totalDlMid3; totalDlMid4];
     totalDlBasal = [totalDlBasal1; totalDlBasal2; totalDlBasal3; totalDlBasal4];
     venusMid = [venusMid1 venusMid2 venusMid3 venusMid4];
     venusBasal = [venusBasal1 venusBasal2 venusBasal3 venusBasal4];
    C = [venusBasal1 venusBasal2 venusBasal3 venusBasal4 ...
        venusMid1 venusMid2 venusMid3 venusMid4];
    D = [totalDlBasal1; totalDlBasal2; totalDlBasal3; totalDlBasal4;...
        totalDlMid1; totalDlMid2; totalDlMid3; totalDlMid4];
    
%     if (mean(D.^2)-mean(D)^2) == 0 
%         Beta = 0; %stops Beta from -> Inf
%         Alpha = mean(C);
%     else
%         Beta = (mean(D.*C)-mean(D)*mean(C))/(mean(D.^2)-mean(D)^2);
%         Alpha = mean(C)-Beta*mean(D);
%     end
Alpha = 0;
Beta = C(1,1)/D(1,1);
    error = sum((C' - Alpha - Beta*D).^2)/length(D);

figure
hold on
plot(Beta*totalDlMid+Alpha,'-r','LineWidth',1.2)
plot(Beta*totalDlBasal+Alpha,'--g','LineWidth',1.2)
plot(venusMid,'-b')
plot(venusBasal,'--b')
xlabel('Time')
ylabel('Amplitude (AU)')
title('Normalized Simulation Data vs. Dorsal-Venus Data')
legend('Ventral (Sim)','Basal (Sim)','Ventral (Venus)','Basal (Venus)','location','Northwest')

end