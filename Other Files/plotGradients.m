function plotGradients(dlNuc,dlCactNuc,Params)

% beta = zeros(5,1);
% phi = beta;
% xi = beta;
% 
% for i = 1:5
%     beta(i) = Params(11+3*(i-1));
%     phi(i) = Params(12+3*(i-1));
%     xi(i) = Params(13+3*(i-1));
% end
% 
% x = linspace(0,1,300);
% Toll10 = beta(1)./(1+(phi(1).*x).^xi(1));
% Toll11 = beta(2)./(1+(phi(2).*x).^xi(2));
% Toll12 = beta(3)./(1+(phi(3).*x).^xi(3));
% Toll13 = beta(4)./(1+(phi(4).*x).^xi(4));
% Toll14 = beta(5)./(1+(phi(5).*x).^xi(5));
% 
figure
% subplot(1,3,1)
% hold on
% plot(x,Toll10,'--','Color',[0 1 1],'LineWidth',1.2)
% plot(x,Toll11,'--','Color',[0 190/255 1],'LineWidth',1.2)
% plot(x,Toll12,'--','Color',[0 150/255 1],'LineWidth',1.2)
% plot(x,Toll13,'--','Color',[0 80/255 1],'LineWidth',1.2)
% plot(x,Toll14,'--','Color',[0 0 1],'LineWidth',1.2)
% title('Toll Signal')
% legend('Toll 10','Toll 11','Toll 12','Toll 13','Toll 14')
    %
    % Greg's data
    %
    [~,N11D,N12D,N13D,N14D] = gregsData;
d11 = N11D(:,end)-min(N11D(:,end));
d11 = (d11/max(d11));
d12 = N12D(:,end)-min(N12D(:,end));
d12 = (d12/max(d12));
d13 = N13D(:,end)-min(N13D(:,end));
d13 = (d13/max(d13));
d14 = N14D(:,end)-min(N14D(:,end));
d14 = (d14/max(d14));



%
% Plotting raw gradients
%
subplot(1,2,1)
hold on
plot(linspace(0,1,length(dlNuc.NC10(:,end))),dlNuc.NC10(:,end)+dlCactNuc.NC10(:,end)-min(dlNuc.NC10(:,end)+dlCactNuc.NC10(:,end)),'-.','Color',[0 1 15/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC11(:,end))),dlNuc.NC11(:,end)+dlCactNuc.NC11(:,end)-min(dlNuc.NC11(:,end)+dlCactNuc.NC11(:,end)),'-.','Color',[0 230/255 14/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC12(:,end))),dlNuc.NC12(:,end)+dlCactNuc.NC12(:,end)-min(dlNuc.NC12(:,end)+dlCactNuc.NC12(:,end)),'-.','Color',[0 200/255 13/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC13(:,end))),dlNuc.NC13(:,end)+dlCactNuc.NC13(:,end)-min(dlNuc.NC13(:,end)+dlCactNuc.NC13(:,end)),'-.','Color',[0 170/255 12/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end)+dlCactNuc.NC14(:,end)-min(dlNuc.NC14(:,end)+dlCactNuc.NC14(:,end)),'-.','Color',[0 130/255 11/255],'LineWidth',1.2)
title('Total dlNuc')
legend('NC10','NC11','NC12','NC13','NC14')
% 
% Plotting normalized gradients
% 
    %
    % Simulation
    %
e10 = dlNuc.NC10(:,end)+dlCactNuc.NC10(:,end)-min(dlNuc.NC10(:,end)+dlCactNuc.NC10(:,end));
e10 = e10/max(e10);
e11 = dlNuc.NC11(:,end)+dlCactNuc.NC11(:,end)-min(dlNuc.NC11(:,end)+dlCactNuc.NC11(:,end));
e11 = e11/max(e11);
e12 = dlNuc.NC12(:,end)+dlCactNuc.NC12(:,end)-min(dlNuc.NC12(:,end)+dlCactNuc.NC12(:,end));
e12 = e12/max(e12);
e13 = dlNuc.NC13(:,end)+dlCactNuc.NC13(:,end)-min(dlNuc.NC13(:,end)+dlCactNuc.NC13(:,end));
e13 = e13/max(e13);
e14 = dlNuc.NC14(:,end)+dlCactNuc.NC14(:,end)-min(dlNuc.NC14(:,end)+dlCactNuc.NC14(:,end));
e14 = e14/max(e14);


    
subplot(1,2,2)
hold on
plot(linspace(0,1,length(dlNuc.NC10(:,end))),e10,'-.','Color',[0 1 15/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC11(:,end))),e11,'-.','Color',[0 230/255 14/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC12(:,end))),e12,'-.','Color',[0 200/255 13/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC13(:,end))),e13,'-.','Color',[0 170/255 12/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC14(:,end))),e14,'-.','Color',[0 130/255 11/255],'LineWidth',1.2)
plot(linspace(0,1,length(N11D(:,end))),d11,'-','Color',[1 130/255 0],'LineWidth',1.2)
plot(linspace(0,1,length(N12D(:,end))),d12,'-','Color',[1 100/255 0],'LineWidth',1.2)
plot(linspace(0,1,length(N13D(:,end))),d13,'-','Color',[1 70/255 0],'LineWidth',1.2)
plot(linspace(0,1,length(N14D(:,end))),d14,'-','Color',[1 40/255 0],'LineWidth',1.2)
title('Total dlNuc Normalized')

% Toll10 = 1./(1+(phi(1).*x).^xi(1));
% Toll11 = 1./(1+(phi(2).*x).^xi(2));
% Toll12 = 1./(1+(phi(3).*x).^xi(3));
% Toll13 = 1./(1+(phi(4).*x).^xi(4));
% Toll14 = 1./(1+(phi(5).*x).^xi(5));

% plot(x,Toll10,'--','Color',[0 1 1],'LineWidth',1.2)
% plot(x,Toll11,'--','Color',[0 190/255 1],'LineWidth',1.2)
% plot(x,Toll12,'--','Color',[0 150/255 1],'LineWidth',1.2)
% plot(x,Toll13,'--','Color',[0 80/255 1],'LineWidth',1.2)
% plot(x,Toll14,'--','Color',[0 0 1],'LineWidth',1.2)
legend('NC10','NC11','NC12','NC13','NC14','Data')

figure
hold on
plot(linspace(0,1,length(dlNuc.NC10(:,end))),e10,'-.','Color',[0 1 15/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC11(:,end))),e11,'-.','Color',[0 230/255 14/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC12(:,end))),e12,'-.','Color',[0 200/255 13/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC13(:,end))),e13,'-.','Color',[0 170/255 12/255],'LineWidth',1.2)
plot(linspace(0,1,length(dlNuc.NC14(:,end))),e14,'-.','Color',[0 130/255 11/255],'LineWidth',1.2)
plot(linspace(0,1,length(N11D(:,end))),d11,'-','Color',[1 130/255 0],'LineWidth',1.2)
plot(linspace(0,1,length(N12D(:,end))),d12,'-','Color',[1 100/255 0],'LineWidth',1.2)
plot(linspace(0,1,length(N13D(:,end))),d13,'-','Color',[1 70/255 0],'LineWidth',1.2)
plot(linspace(0,1,length(N14D(:,end))),d14,'-','Color',[1 40/255 0],'LineWidth',1.2)

plot(linspace(0,-1,length(dlNuc.NC10(:,end))),e10,'-.','Color',[0 1 15/255],'LineWidth',1.2)
plot(linspace(0,-1,length(dlNuc.NC11(:,end))),e11,'-.','Color',[0 230/255 14/255],'LineWidth',1.2)
plot(linspace(0,-1,length(dlNuc.NC12(:,end))),e12,'-.','Color',[0 200/255 13/255],'LineWidth',1.2)
plot(linspace(0,-1,length(dlNuc.NC13(:,end))),e13,'-.','Color',[0 170/255 12/255],'LineWidth',1.2)
plot(linspace(0,-1,length(dlNuc.NC14(:,end))),e14,'-.','Color',[0 130/255 11/255],'LineWidth',1.2)
plot(linspace(0,-1,length(N11D(:,end))),d11,'-','Color',[1 130/255 0],'LineWidth',1.2)
plot(linspace(0,-1,length(N12D(:,end))),d12,'-','Color',[1 100/255 0],'LineWidth',1.2)
plot(linspace(0,-1,length(N13D(:,end))),d13,'-','Color',[1 70/255 0],'LineWidth',1.2)
plot(linspace(0,-1,length(N14D(:,end))),d14,'-','Color',[1 40/255 0],'LineWidth',1.2)

title('Total dlNuc Normalized')


figure
hold on
plot(linspace(0,1,length(dlNuc.NC10(:,end))),dlNuc.NC10(:,end)+dlCactNuc.NC10(:,end),'-.','Color','black','LineWidth',2)
plot(linspace(0,1,length(dlNuc.NC11(:,end))),dlNuc.NC11(:,end)+dlCactNuc.NC11(:,end),'Color',[35 68 252]/255,'LineWidth',2)
plot(linspace(0,1,length(dlNuc.NC12(:,end))),dlNuc.NC12(:,end)+dlCactNuc.NC12(:,end),'Color',[47 201 224]/255,'LineWidth',2)
plot(linspace(0,1,length(dlNuc.NC13(:,end))),dlNuc.NC13(:,end)+dlCactNuc.NC13(:,end),'Color',[209 194 27]/255,'LineWidth',2)
plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end)+dlCactNuc.NC14(:,end),'Color',[245 39 49]/255,'LineWidth',2)

plot(linspace(0,-1,length(dlNuc.NC10(:,end))),dlNuc.NC10(:,end)+dlCactNuc.NC10(:,end),'-.','Color','black','LineWidth',2)
plot(linspace(0,-1,length(dlNuc.NC11(:,end))),dlNuc.NC11(:,end)+dlCactNuc.NC11(:,end),'Color',[35 68 252]/255,'LineWidth',2)
plot(linspace(0,-1,length(dlNuc.NC12(:,end))),dlNuc.NC12(:,end)+dlCactNuc.NC12(:,end),'Color',[47 201 224]/255,'LineWidth',2)
plot(linspace(0,-1,length(dlNuc.NC13(:,end))),dlNuc.NC13(:,end)+dlCactNuc.NC13(:,end),'Color',[209 194 27]/255,'LineWidth',2)
plot(linspace(0,-1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end)+dlCactNuc.NC14(:,end),'Color',[245 39 49]/255,'LineWidth',2)
ylim([0 3.5])
title('Total dlNuc')
legend('NC10','NC11','NC12','NC13','NC14')
end