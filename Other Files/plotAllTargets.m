function plotAllTargets(protein,time)

dlNuc = protein{1};
dlCactNuc = protein{3};

snaCyt = protein{7};
sogCyt = protein{8};
dppCyt = protein{9};
vndCyt = protein{10};

% 3D Profiles
figure
subplot(2,2,1)
surf(time.T10',time.X10',snaCyt.NC10)
hold on
surf(time.T10m',time.X10m',snaCyt.NC10m)
surf(time.T11',time.X11',snaCyt.NC11)
surf(time.T11m',time.X11m',snaCyt.NC11m)
surf(time.T12',time.X12',snaCyt.NC12)
surf(time.T12m',time.X12m',snaCyt.NC12m)
surf(time.T13',time.X13',snaCyt.NC13)
surf(time.T13m',time.X13m',snaCyt.NC13m)
surf(time.T14',time.X14',snaCyt.NC14)
ylabel('space')
xlabel('time')
title('Snail')
shading flat


subplot(2,2,2)
surf(time.T10',time.X10',sogCyt.NC10)
hold on
surf(time.T10m',time.X10m',sogCyt.NC10m)
surf(time.T11',time.X11',sogCyt.NC11)
surf(time.T11m',time.X11m',sogCyt.NC11m)
surf(time.T12',time.X12',sogCyt.NC12)
surf(time.T12m',time.X12m',sogCyt.NC12m)
surf(time.T13',time.X13',sogCyt.NC13)
surf(time.T13m',time.X13m',sogCyt.NC13m)
surf(time.T14',time.X14',sogCyt.NC14)
ylabel('space')
xlabel('time')
title('Sog')
shading flat


subplot(2,2,3)
surf(time.T10',time.X10',dppCyt.NC10)
hold on
surf(time.T10m',time.X10m',dppCyt.NC10m)
surf(time.T11',time.X11',dppCyt.NC11)
surf(time.T11m',time.X11m',dppCyt.NC11m)
surf(time.T12',time.X12',dppCyt.NC12)
surf(time.T12m',time.X12m',dppCyt.NC12m)
surf(time.T13',time.X13',dppCyt.NC13)
surf(time.T13m',time.X13m',dppCyt.NC13m)
surf(time.T14',time.X14',dppCyt.NC14)
ylabel('space')
xlabel('time')
title('Dpp')
shading flat


subplot(2,2,4)
surf(time.T10',time.X10',vndCyt.NC10)
hold on
surf(time.T10m',time.X10m',vndCyt.NC10m)
surf(time.T11',time.X11',vndCyt.NC11)
surf(time.T11m',time.X11m',vndCyt.NC11m)
surf(time.T12',time.X12',vndCyt.NC12)
surf(time.T12m',time.X12m',vndCyt.NC12m)
surf(time.T13',time.X13',vndCyt.NC13)
surf(time.T13m',time.X13m',vndCyt.NC13m)
surf(time.T14',time.X14',vndCyt.NC14)
ylabel('space')
xlabel('time')
title('Vnd')
shading flat


x = linspace(0,1,51);

dlNorm = (dlNuc.NC14(:,end)+dlCactNuc.NC14(:,end)-min(dlNuc.NC14(:,end)+...
    dlCactNuc.NC14(:,end)));
dlNorm = dlNorm/max(dlNorm);

% Normalized 2D profiles
figure
plot(x,dlNorm,'--','LineWidth',3,'Color','black')
hold on
plot(x,sogCyt.NC14(:,end)/max(sogCyt.NC14(:,end)),'LineWidth',2.2,'Color','g')
plot(x,snaCyt.NC14(:,end)/max(snaCyt.NC14(:,end)),'LineWidth',2.2,'Color','r')
plot(x,dppCyt.NC14(:,end)/max(dppCyt.NC14(:,end)),'LineWidth',2.2,'Color','y')
plot(x,vndCyt.NC14(:,end)/max(vndCyt.NC14(:,end)),'LineWidth',2.2,'Color','b')
title('Normalized mRNA Expression (Late NC14)')
legend('Dorsal','sog','sna','dpp','vnd')
xlabel('DV Axis Location')
ylabel('Normalized Amplitude (AU)')

end