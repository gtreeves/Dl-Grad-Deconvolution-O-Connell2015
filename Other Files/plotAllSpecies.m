function plotAllSpecies(Model,protein,time)

% dlNuc = protein{1};
% dlCyt = protein{2};
% dlCactNuc = protein{3};
% dlCactCyt = protein{4};
% cactNuc = protein{5};
% cactCyt = protein{6};

switch Model
    case 'shuttling'
        names = {'dlNuc','dlCyt','dlCactCyt','cactCyt'};
        figure
        hold on
        for i = 1:4
            species=protein{i};
            
            subplot(2,2,i)
            surf(time.T10',time.X10',species.NC10)
            hold on
            surf(time.T10m',time.X10m',species.NC10m)
            surf(time.T11',time.X11',species.NC11)
            surf(time.T11m',time.X11m',species.NC11m)
            surf(time.T12',time.X12',species.NC12)
            surf(time.T12m',time.X12m',species.NC12m)
            surf(time.T13',time.X13',species.NC13)
            surf(time.T13m',time.X13m',species.NC13m)
            surf(time.T14',time.X14',species.NC14)
            surf(time.T10',-time.X10',species.NC10)
            surf(time.T10m',-time.X10m',species.NC10m)
            surf(time.T11',-time.X11',species.NC11)
            surf(time.T11m',-time.X11m',species.NC11m)
            surf(time.T12',-time.X12',species.NC12)
            surf(time.T12m',-time.X12m',species.NC12m)
            surf(time.T13',-time.X13',species.NC13)
            surf(time.T13m',-time.X13m',species.NC13m)
            surf(time.T14',-time.X14',species.NC14)
            ylabel('DV Coordinate ')
            xlabel('time (min) ')
            title(names{i})
            shading flat
            set(gcf,'Position',[1,1,350,350])
            
%             figure
%             hold on
%             nc11norm = nc11-(min(nc11));
%             nc11norm = nc11norm/max(nc11norm);
%             plot(linspace(0,1,length(nc11norm)),nc11norm)
%             nc12norm = nc11-(min(nc12));
%             nc12norm = nc12norm/max(nc12norm);
%             plot(linspace(0,1,length(nc12norm)),nc12norm)
%             nc13norm = nc11-(min(nc13));
%             nc13norm = nc13norm/max(nc13norm);
%             plot(linspace(0,1,length(nc13norm)),nc13norm)
%             nc14norm = nc14-(min(nc14));
%             nc14norm = nc14norm/max(nc14norm);
%             plot(linspace(0,1,length(nc14norm)),nc14norm)
        end
        
    case 'tollSaturation'
        pIndex = [1 2 4 6 7 8]; % which cells in 'protein' do we need?
        names = {'dlNuc','dlCyt','dlCactCyt','cactCyt','tollActive','TDC-complex'};
        figure
        for i = 1:6
            species=protein{pIndex(i)};
            
            subplot(2,3,i)
            surf(time.T10',time.X10',species.NC10)
            hold on
            surf(time.T10m',time.X10m',species.NC10m)
            surf(time.T11',time.X11',species.NC11)
            surf(time.T11m',time.X11m',species.NC11m)
            surf(time.T12',time.X12',species.NC12)
            surf(time.T12m',time.X12m',species.NC12m)
            surf(time.T13',time.X13',species.NC13)
            surf(time.T13m',time.X13m',species.NC13m)
            surf(time.T14',time.X14',species.NC14)
            ylabel('space')
            xlabel('time')
            title(names{i})
            shading flat
        end
        
        
    case 'dynamicDorsal'
        dlNuc = protein{1};
        dlCyt = protein{2};
        dlCactNuc = protein{3};
        dlCactCyt = protein{4};
        cactNuc = protein{5};
        cactCyt = protein{6};
        names = {'dlNuc','dlCyt','dlCactNuc','dlCactCyt','cactNuc','cactCyt'};
        figure
        for i = 1:6
            species=protein{i};
            subplot(2,3,i)
            surf(time.T10',time.X10',species.NC10)
            hold on
            surf(time.T10m',time.X10m',species.NC10m)
            surf(time.T11',time.X11',species.NC11)
            surf(time.T11m',time.X11m',species.NC11m)
            surf(time.T12',time.X12',species.NC12)
            surf(time.T12m',time.X12m',species.NC12m)
            surf(time.T13',time.X13',species.NC13)
            surf(time.T13m',time.X13m',species.NC13m)
            surf(time.T14',time.X14',species.NC14)
            ylabel('space')
            xlabel('time')
            title(names{i})
            shading flat
        end
        
        %% Sum of free and complexed dorsal
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
        
        %% Free nuclear dorsal
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
        title('Free nuclear dl')
        shading flat
        
        
        %% Complexed nuclear dorsal
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
        title('Complexed nuclear dl')
        shading flat
        
        %% Distinguishing free and complexed dorsal at NC14
        figure
        plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end),...
            'b','LineWidth',3)
        hold on
        plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end)...
            +dlCactNuc.NC14(:,end),'r','LineWidth',3)
        title('Distinguishing Free and Complexed Dorsal','FontSize',16)
        xlabel('DV Coordinate','FontSize',16)
        ylabel('Amplitude (AU) ','FontSize',16)
        legend('Free dl','Free+Complexed')
        
        %% Distinguishing free and complexed dorsal at all NCs
        figure
        subplot(1,2,1)
        plot(linspace(0,1,length(dlNuc.NC11(:,end))),dlNuc.NC11(:,end),...
            'b','LineWidth',1)
        hold on
        plot(linspace(0,1,length(dlNuc.NC12(:,end))),dlNuc.NC12(:,end),...
            'b','LineWidth',1.5)
        plot(linspace(0,1,length(dlNuc.NC13(:,end))),dlNuc.NC13(:,end),...
            'b','LineWidth',2)
        plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end),...
            'b','LineWidth',3)
        title('Free Dorsal')
        
        subplot(1,2,2)
        plot(linspace(0,1,length(dlNuc.NC11(:,end))),dlNuc.NC11(:,end)...
            +dlCactNuc.NC11(:,end),'r','LineWidth',1)
        hold on
        plot(linspace(0,1,length(dlNuc.NC12(:,end))),dlNuc.NC12(:,end)...
            +dlCactNuc.NC12(:,end),'r','LineWidth',1.5)
        plot(linspace(0,1,length(dlNuc.NC13(:,end))),dlNuc.NC13(:,end)...
            +dlCactNuc.NC13(:,end),'r','LineWidth',2)
        plot(linspace(0,1,length(dlNuc.NC14(:,end))),dlNuc.NC14(:,end)...
            +dlCactNuc.NC14(:,end),'r','LineWidth',3)
        title('Free + Complexed Dorsal')
        
end
end