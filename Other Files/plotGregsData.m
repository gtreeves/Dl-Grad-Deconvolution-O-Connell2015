function [varargout]=plotGregsData

%% Greg's data
load('hybrid_embryo.mat')

A = data.A;
B = data.B;
t = data.t;
SIG = data.Sig;

%% NC 11
N11T = t(1:16);
N11A = A(1:16);
N11B = B(1:16);
sig = 0.15;
m = -0.07;
M = 25;
x = linspace(0,1,M);
N11D = repmat(N11A,1,M).*exp(-repmat(x,length(N11A),1).^2/2/sig.^2) + ...
    repmat(N11B,1,M) + m*N11A*x;
N11D = fliplr(rot90(N11D,3));

%% NC 12
N12T = t(33:60);
N12A = A(33:60);
N12B = B(33:60);

x = linspace(0,1,M);
N12D = repmat(N12A,1,M).*exp(-repmat(x,length(N12A),1).^2/2/sig.^2) + ...
    repmat(N12B,1,M) + m*N12A*x;
N12D = fliplr(rot90(N12D,3));

%% NC 13
N13T = t(76:126);
N13A = A(76:126);
N13B = B(76:126);
x = linspace(0,1,M);
N13D = repmat(N13A,1,M).*exp(-repmat(x,length(N13A),1).^2/2/sig.^2) + ...
    repmat(N13B,1,M) + m*N13A*x;
N13D = fliplr(rot90(N13D,3));

%% NC 14
N14T = t(147:327);
N14A = A(147:327);
N14B = B(147:327);
x = linspace(0,1,M);
N14D = repmat(N14A,1,M).*exp(-repmat(x,length(N14A),1).^2/2/sig.^2) + ...
    repmat(N14B,1,M) + m*N14A*x;
N14D = fliplr(rot90(N14D,3));

% figure
% hold on
% surf(N11T,x,N11D)
% surf(N11T,-x,N11D)
% surf(N12T,x,N12D)
% surf(N12T,-x,N12D)
% surf(N13T,x,N13D)
% surf(N13T,-x,N13D)
% surf(N14T,x,N14D)
% surf(N14T,-x,N14D)
% shading flat


allData = repmat(A,1,M).*exp(-repmat(x,length(A),1).^2/2/sig.^2) + ...
    repmat(B,1,M) + m*A*x;
A2 = data.A2;
B2 = data.B2;
Sig2 = data.Sig2;

totalDorsalAllData = repmat(A2,1,M).*exp(-repmat(x,length(A2),1).^2./2./repmat(Sig2,1,25).^2) + ...
    repmat(B2,1,M) + m*A2*x;
% figure
% hold on
% surf(x,t,allData)
% surf(-x,t,allData)
% shading flat
% view(78,20)

if nargout == 0
% figure
% hold on
% surf(x,t(1:375)+7.695,allData(1:375,:)/max(max(allData(1:375,:))))
% surf(-x,t(1:375)+7.695,allData(1:375,:)/max(max(allData(1:375,:))))
% shading flat
% view(42,42)
% title('Nuclear dl fluorescence (normalized)')

figure
hold on
surf(x,t(1:375)+7.695,allData(1:375,:))
surf(-x,t(1:375)+7.695,allData(1:375,:))
plot3(zeros(375,1),t(1:375)+7.695,allData(1:375,1),'Color','blue','LineWidth',3)
plot3(-1*ones(375,1),t(1:375)+7.695,allData(1:375,end),'Color','red','LineWidth',3)
plot3(ones(375,1),t(1:375)+7.695,allData(1:375,end),'Color','red','LineWidth',3)

shading flat
view(42,42)

% title('Nuclear dl fluorescence','FontName','Arial','FontSize',12)
xlabel('DV Coordinate','FontName','Arial','FontSize',12)
ylabel('time (min)','FontName','Arial','FontSize',12)
zlabel('Intensity (AU)','FontName','Arial','FontSize',12)
set(gca,'FontName','Arial','FontSize',12)
grid on
% 
% figure
% hold on
% plot(x,N11D(:,end))
% plot(-x,N11D(:,end))
% plot(x,N12D(:,end))
% plot(-x,N12D(:,end))
% plot(x,N13D(:,end))
% plot(-x,N13D(:,end))
% plot(x,N14D(:,end))
% plot(-x,N14D(:,end))
% shading flat

% figure
% hold on
% surf(x,t(1:375),totalDorsalAllData(1:375,:)/max(max(totalDorsalAllData(1:375,:))))
% surf(-x,t(1:375),totalDorsalAllData(1:375,:)/max(max(totalDorsalAllData(1:375,:))))
% shading flat
% view(42,42)
% title('Total dl fluorescence (normalized)')
% 
% figure
% hold on
% surf(x,t(1:375),allData(1:375,:)-totalDorsalAllData(1:375,:))
% surf(-x,t(1:375),allData(1:375,:)-totalDorsalAllData(1:375,:))
% shading flat
% view(42,42)
% title('Cytoplasmic dl fluorescence (normalized)?')
% plot3(zeros(375,1),t(1:375),allData(1:375,1),'LineWidth',3,'Color','blue')
% plot3(ones(375,1),t(1:375),allData(1:375,end),'LineWidth',3,'Color','red')
% plot3(-ones(375,1),t(1:375),allData(1:375,end),'LineWidth',3,'Color','red')
end




dat = {allData,totalDorsalAllData,x,t,SIG,N11D,N12D,N13D,N14D,N11T,N12T,N13T,N14T};
% assign varargout
varargout = cell(nargout);
for i = 1:max(nargout,1)
    varargout{i} = dat{i};
end

end



