function plotDlCactSpecies(protein,time)

dlNuc = protein{1};
dlCyt = protein{2};
dlCactNuc = protein{3};
dlCactCyt = protein{4};
cactNuc = protein{5};
cactCyt = protein{6};

figure
subplot(2,3,1)
surf(time.T10',time.X10',dlNuc.NC10)
hold on
surf(time.T11',time.X11',dlNuc.NC11)
surf(time.T12',time.X12',dlNuc.NC12)
surf(time.T13',time.X13',dlNuc.NC13)
surf(time.T14',time.X14',dlNuc.NC14)
ylabel('space')
xlabel('time')
title('dlNuc')
shading flat

subplot(2,3,4)
surf(time.T10',time.X10',dlCyt.NC10)
hold on
surf(time.T10m',time.X10m',dlCyt.NC10m)
surf(time.T11',time.X11',dlCyt.NC11)
surf(time.T11m',time.X11m',dlCyt.NC11m)
surf(time.T12',time.X12',dlCyt.NC12)
surf(time.T12m',time.X12m',dlCyt.NC12m)
surf(time.T13',time.X13',dlCyt.NC13)
surf(time.T13m',time.X13m',dlCyt.NC13m)
surf(time.T14',time.X14',dlCyt.NC14)
ylabel('space')
xlabel('time')
title('dlCyt')
shading flat

subplot(2,3,2)
surf(time.T10',time.X10',cactNuc.NC10)
hold on
surf(time.T11',time.X11',cactNuc.NC11)
surf(time.T12',time.X12',cactNuc.NC12)
surf(time.T13',time.X13',cactNuc.NC13)
surf(time.T14',time.X14',cactNuc.NC14)
ylabel('space')
xlabel('time')
title('cactNuc')
shading flat

subplot(2,3,5)
surf(time.T10',time.X10',cactCyt.NC10)
hold on
surf(time.T10m',time.X10m',cactCyt.NC10m)
surf(time.T11',time.X11',cactCyt.NC11)
surf(time.T11m',time.X11m',cactCyt.NC11m)
surf(time.T12',time.X12',cactCyt.NC12)
surf(time.T12m',time.X12m',cactCyt.NC12m)
surf(time.T13',time.X13',cactCyt.NC13)
surf(time.T13m',time.X13m',cactCyt.NC13m)
surf(time.T14',time.X14',cactCyt.NC14)
ylabel('space')
xlabel('time')
title('cactCyt')
shading flat

subplot(2,3,3)
surf(time.T10',time.X10',dlCactNuc.NC10)
hold on
surf(time.T11',time.X11',dlCactNuc.NC11)
surf(time.T12',time.X12',dlCactNuc.NC12)
surf(time.T13',time.X13',dlCactNuc.NC13)
surf(time.T14',time.X14',dlCactNuc.NC14)
ylabel('space')
xlabel('time')
title('dlCactNuc')
shading flat

subplot(2,3,6)
surf(time.T10',time.X10',dlCactCyt.NC10)
hold on
surf(time.T10m',time.X10m',dlCactCyt.NC10m)
surf(time.T11',time.X11',dlCactCyt.NC11)
surf(time.T11m',time.X11m',dlCactCyt.NC11m)
surf(time.T12',time.X12',dlCactCyt.NC12)
surf(time.T12m',time.X12m',dlCactCyt.NC12m)
surf(time.T13',time.X13',dlCactCyt.NC13)
surf(time.T13m',time.X13m',dlCactCyt.NC13m)
surf(time.T14',time.X14',dlCactCyt.NC14)
ylabel('space')
xlabel('time')
title('dlCactCyt')
shading flat

figure
surf(time.T10',time.X10',dlNuc.NC10+dlCactNuc.NC10)
hold on
surf(time.T11',time.X11',dlNuc.NC11+dlCactNuc.NC11)
surf(time.T12',time.X12',dlNuc.NC12+dlCactNuc.NC12)
surf(time.T13',time.X13',dlNuc.NC13+dlCactNuc.NC13)
surf(time.T14',time.X14',dlNuc.NC14+dlCactNuc.NC14)

surf(time.T10',-time.X10',dlNuc.NC10+dlCactNuc.NC10)
surf(time.T11',-time.X11',dlNuc.NC11+dlCactNuc.NC11)
surf(time.T12',-time.X12',dlNuc.NC12+dlCactNuc.NC12)
surf(time.T13',-time.X13',dlNuc.NC13+dlCactNuc.NC13)
surf(time.T14',-time.X14',dlNuc.NC14+dlCactNuc.NC14)

ylabel('space')
xlabel('time')
title('Total nuclear dl (free & complexed)')
shading flat

figure
surf(time.T10',time.X10',dlNuc.NC10)
hold on
surf(time.T11',time.X11',dlNuc.NC11)
surf(time.T12',time.X12',dlNuc.NC12)
surf(time.T13',time.X13',dlNuc.NC13)
surf(time.T14',time.X14',dlNuc.NC14)

surf(time.T10',-time.X10',dlNuc.NC10)
surf(time.T11',-time.X11',dlNuc.NC11)
surf(time.T12',-time.X12',dlNuc.NC12)
surf(time.T13',-time.X13',dlNuc.NC13)
surf(time.T14',-time.X14',dlNuc.NC14)

ylabel('space')
xlabel('time')
title('Total nuclear dl (free & complexed)')
shading flat

figure
surf(time.T10',time.X10',dlCactNuc.NC10)
hold on
surf(time.T11',time.X11',dlCactNuc.NC11)
surf(time.T12',time.X12',dlCactNuc.NC12)
surf(time.T13',time.X13',dlCactNuc.NC13)
surf(time.T14',time.X14',dlCactNuc.NC14)

surf(time.T10',-time.X10',dlCactNuc.NC10)
surf(time.T11',-time.X11',dlCactNuc.NC11)
surf(time.T12',-time.X12',dlCactNuc.NC12)
surf(time.T13',-time.X13',dlCactNuc.NC13)
surf(time.T14',-time.X14',dlCactNuc.NC14)

ylabel('space')
xlabel('time')
title('Total nuclear dl (free & complexed)')
shading flat

figure
plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end),'b','LineWidth',3)
hold on
plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end)+dlCactNuc.NC14(:,end),'r','LineWidth',3)
title('Distinguishing Free and Complexed Dorsal','FontSize',16)
xlabel('DV Coordinate','FontSize',16)
ylabel('Amplitude (AU) ','FontSize',16)

figure
subplot(1,2,1)
plot(linspace(0,1,length(dlNuc.NC11(:,end))),dlNuc.NC11(:,end),'b','LineWidth',3)
hold on
plot(linspace(0,1,length(dlNuc.NC12(:,end))),dlNuc.NC12(:,end),'b','LineWidth',3)
plot(linspace(0,1,length(dlNuc.NC13(:,end))),dlNuc.NC13(:,end),'b','LineWidth',3)
plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end),'b','LineWidth',3)

subplot(1,2,2)
plot(linspace(0,1,length(dlNuc.NC11(:,end))),dlNuc.NC11(:,end)+dlCactNuc.NC11(:,end),'r','LineWidth',3)
hold on
plot(linspace(0,1,length(dlNuc.NC12(:,end))),dlNuc.NC12(:,end)+dlCactNuc.NC12(:,end),'r','LineWidth',3)
plot(linspace(0,1,length(dlNuc.NC13(:,end))),dlNuc.NC13(:,end)+dlCactNuc.NC13(:,end),'r','LineWidth',3)
plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end)+dlCactNuc.NC14(:,end),'r','LineWidth',3)


end