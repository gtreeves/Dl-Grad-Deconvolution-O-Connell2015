function plotData(protein,time,plot_type,varargin)

ncs = {'NC10','NC10m','NC11','NC11m','NC12','NC12m','NC13','NC13m','NC14'};
ts = {'T10','T10m','T11','T11m','T12','T12m','T13','T13m','T14'};
xs = {'X10','X10m','X11','X11m','X12','X12m','X13','X13m','X14'};

Beta = 1;
figure_type = 'gui';

for i = 1:2:length(varargin)
    argName = lower(varargin{i});
    switch argName
        case 'beta'
            Beta = varargin{i+1};
        case 'figure_type'
            figure_type = varargin{i+1};
    end
end


if strcmp(figure_type,'new2Dfigure')
    figure
    hold on
    switch plot_type
        case 'nc14 gene expression'
            plotGeneExpression2D;
        case 'nuclearDorsal'
            plotNuclearDorsal2D;
        otherwise        
            plotProtein2D;
    end
elseif strcmp(figure_type,'new3Dfigure')
    figure
    hold on
    switch plot_type % for generating plots in a new window
        case {'nc14 gene expression'}
            plotGeneExpression3D;
        case 'nuclearDorsal'
            plotNuclearDorsal3D;
        otherwise
            plotProtein3D;
    end
else %for plotting in the gui
    switch plot_type
        case 'nc14 gene expression'
            plotGeneExpression2D;
            
        otherwise
            plotProtein3D;
            
    end
end


    function plotGeneExpression2D
        hold off
            plot(linspace(0,1,51),protein.sna.NC14(:,end),'Color',[0.6 0 0],'LineWidth',2)
            hold on
            plot(linspace(0,1,51),protein.sog.NC14(:,end),'Color',[31 184 50]/255,'LineWidth',2)
            plot(linspace(0,1,51),protein.vnd.NC14(:,end),'Color',[38 65 240]/255,'LineWidth',2)
            plot(linspace(0,1,51),protein.dpp.NC14(:,end),'Color',[240 213  38]/255,'LineWidth',2)
            xlabel('DV Coordinate ')
            ylabel('expression ')
            set(gcf,'Name','mRNA_fig')
    end

    function plotNuclearDorsal2D
        [allData,~,~,t,sigma] = plotGregsData2;
%             if width
%                 [Sig,A,B,M,SIG] = gradientWidth(protein.dlNuc,protein.dlCactNuc);
%             end
            dataAmp = allData(1:375,1)-allData(1:375,end);
            
            rectangle('Position',[7.695,0,3.605,1],'FaceColor',[230 230 230]/255,'LineStyle','None')
            rectangle('Position',[15.63,0,6.01,1],'FaceColor',[230 230 230]/255,'LineStyle','None')
            rectangle('Position',[25.97,0,11.79,1],'FaceColor',[230 230 230]/255,'LineStyle','None')
            rectangle('Position',[43.3,0,54.4,1],'FaceColor',[230 230 230]/255,'LineStyle','None')
            
            ncs = {'NC11','NC12','NC13','NC14'};
            ts = {'T11','T12','T13','T14'};
            
            %amplitude
            plot(t(1:375)+7.695,dataAmp/max(dataAmp),'LineWidth',2,...
                'Color',[98 12 232]/255,'LineStyle','-.')
            
            %basal
            plot(t(1:375)+7.695,allData(1:375,end)/max(dataAmp),'LineWidth',2,...
                'Color',[0 255 255]/255,'LineStyle','-.')
            
            %sigma
            plot(t(1:375)+7.695,sigma(1:375),'LineWidth',2,...
                'Color',[255 0 0]/255,'LineStyle','-.')
            
            for i = 1:4
                plot(time.(ts{i})(:,1),Beta*(protein.dlNuc.(ncs{i})(1,:)+protein.dlCactNuc.(ncs{i})(1,:)...
                    -(protein.dlNuc.(ncs{i})(end,:)+protein.dlCactNuc.(ncs{i})(end,:)))/max(dataAmp)...
                    ,'Color',[209 0 255]/255,'LineWidth',3)
                plot(time.(ts{i})(:,1),Beta*(protein.dlNuc.(ncs{i})(end,:)+protein.dlCactNuc.(ncs{i})(end,:))...
                    /max(dataAmp),'Color',[12 111 232]/255,'LineWidth',3)
%                 if width
%                     plot(time.(ts{i})(:,1),Sig.(ncs{i}),'Color','red','LineWidth',2)
%                 end
            end
            
            plot([7.695 100],[0.15 0.15])
            title('Amplitude & Basal Levels','FontSize',12,'FontName','Arial')
            ylabel('Normalized Value' ,'FontSize',12,'FontName','Arial')
            xlabel('time (min) ','FontSize',12,'FontName','Arial')
            ylim([0 1])
            xlim([7 98])
            set(gcf,'Name','Model_results_in_time')
            set(gca,'FontSize',12,'FontName','Arial')
            box on
            
            pause(0.1)
            % from fitgauss
            figure
            hold on
            rectangle('Position',[7.695,0,3.605,1],'FaceColor',[230 230 230]/255,'LineStyle','None')
            rectangle('Position',[15.63,0,6.01,1],'FaceColor',[230 230 230]/255,'LineStyle','None')
            rectangle('Position',[25.97,0,11.79,1],'FaceColor',[230 230 230]/255,'LineStyle','None')
            rectangle('Position',[43.3,0,54.4,1],'FaceColor',[230 230 230]/255,'LineStyle','None')
            
            
            %amplitude
            plot(t(1:375)+7.695,dataAmp/max(dataAmp),'LineWidth',2,...
                'Color',[98 12 232]/255,'LineStyle','-.')
            
            %basal
            plot(t(1:375)+7.695,allData(1:375,end)/max(dataAmp),'LineWidth',2,...
                'Color',[0 255 255]/255,'LineStyle','-.')
            
            for i = 1:4
                plot(time.(ts{i})(:,1),(protein.dlNuc.(ncs{i})(1,:)...
                    +protein.dlCactNuc.(ncs{i})(1,:)-(protein.dlNuc.(ncs{i})(end,:)...
                    +protein.dlCactNuc.(ncs{i})(end,:)))*Beta/max(dataAmp),'Color',...
                    [209 0 255]/255,'LineWidth',3)
                plot(time.(ts{i})(:,1),(protein.dlNuc.(ncs{i})(end,:)...
                    +protein.dlCactNuc.(ncs{i})(end,:))*Beta/max(dataAmp)...
                    ,'Color',[12 111 232]/255,'LineWidth',3)
%                 if width
%                     plot(time.(ts{i})(:,1),SIG.(ncs{i}),'Color','red','LineWidth',2)
%                     %             plot(time.(ts{i})(:,1),M.(ncs{i}),'Color','green','LineWidth',2)
%                 end
            end
            plot([7.695 100],[0.15 0.15])
            title('Amplitude & Basal Levels (Gaussian)','FontSize',12,'FontName','Arial')
            ylabel('Normalized Value' ,'FontSize',12,'FontName','Arial')
            xlabel('time (min) ','FontSize',12,'FontName','Arial')
            
            xlim([7 98])
            set(gcf,'Name','Model_results_in_time_Gaussian')
            set(gca,'FontSize',12,'FontName','Arial')
            box on
            % text(5,500,strcat('Alpha = ',num2str(Alpha),'; Beta = ',num2str(Beta)))
            
            
            figure
            plot(linspace(0,1,25),allData(300,:),'LineWidth',1,'Color',[209 0 255]/255)
            hold on
            plot(linspace(0,-1,25),allData(300,:),'LineWidth',1,'Color',[209 0 255]/255)
            plot(linspace(0,1,25),allData(127,:),'LineWidth',1,'Color',[98 12 232]/255)
            plot(linspace(0,-1,25),allData(127,:),'LineWidth',1,'Color',[98 12 232]/255)
            plot(linspace(0,1,25),allData(62,:),'LineWidth',1,'Color',[12 111 232]/255)
            plot(linspace(0,-1,25),allData(62,:),'LineWidth',1,'Color',[12 111 232]/255)
            plot(linspace(0,1,25),allData(16,:),'LineWidth',1,'Color',[0 255 255]/255)
            plot(linspace(0,-1,25),allData(16,:),'LineWidth',1,'Color',[0 255 255]/255)
            
            plot(linspace(0,1,51),Beta*(protein.dlNuc.NC14(:,end)+protein.dlCactNuc.NC14(:,end)),'LineWidth',2,'Color',[209 0 255]/255)
            plot(linspace(0,-1,51),Beta*(protein.dlNuc.NC14(:,end)+protein.dlCactNuc.NC14(:,end)),'LineWidth',2,'Color',[209 0 255]/255)
            plot(linspace(0,1,36),Beta*(protein.dlNuc.NC13(:,end)+protein.dlCactNuc.NC13(:,end)),'LineWidth',2,'Color',[98 12 232]/255)
            plot(linspace(0,-1,36),Beta*(protein.dlNuc.NC13(:,end)+protein.dlCactNuc.NC13(:,end)),'LineWidth',2,'Color',[98 12 232]/255)
            plot(linspace(0,1,26),Beta*(protein.dlNuc.NC12(:,end)+protein.dlCactNuc.NC12(:,end)),'LineWidth',2,'Color',[12 111 232]/255)
            plot(linspace(0,-1,26),Beta*(protein.dlNuc.NC12(:,end)+protein.dlCactNuc.NC12(:,end)),'LineWidth',2,'Color',[12 111 232]/255)
            plot(linspace(0,1,19),Beta*(protein.dlNuc.NC11(:,end)+protein.dlCactNuc.NC11(:,end)),'LineWidth',2,'Color',[0 255 255]/255)
            plot(linspace(0,-1,19),Beta*(protein.dlNuc.NC11(:,end)+protein.dlCactNuc.NC11(:,end)),'LineWidth',2,'Color',[0 255 255]/255)
            title('Gradient Comparison','FontSize',12,'FontName','Arial')
            xlabel('DV Coordinate ','FontSize',12,'FontName','Arial')
            ylabel('Intensity (AU)','FontSize',12,'FontName','Arial')
            ylim([0 6000])
            set(gcf,'Name','Gradient_comparison_in_space')
            set(gca,'FontSize',12,'FontName','Arial')
            box on
            
            n11=(protein.dlNuc.NC11(:,end)+protein.dlCactNuc.NC11(:,end));
            n11 = (n11-min(n11))/max(n11-min(n11))';
            n12=(protein.dlNuc.NC12(:,end)+protein.dlCactNuc.NC12(:,end));
            n12 = (n12-min(n12))/max(n12-min(n12))';
            n13=(protein.dlNuc.NC13(:,end)+protein.dlCactNuc.NC13(:,end));
            n13 = (n13-min(n13))/max(n13-min(n13))';
            n14=(protein.dlNuc.NC14(:,end)+protein.dlCactNuc.NC14(:,end));
            n14 = (n14-min(n14))/max(n14-min(n14))';
            
            d11=(allData(16,:)-min(allData(16,:)))/max(allData(16,:)-min(allData(16,:)));
            d12=(allData(62,:)-min(allData(62,:)))/max(allData(62,:)-min(allData(62,:)));
            d13=(allData(127,:)-min(allData(127,:)))/max(allData(127,:)-min(allData(127,:)));
            d14=(allData(300,:)-min(allData(300,:)))/max(allData(300,:)-min(allData(300,:)));
            
            
            figure
            plot([linspace(-1,0,25) linspace(0,1,25)],[fliplr(d14) d14],'LineWidth',1,'Color',[209 0 255]/255)
            hold on
            plot([linspace(-1,0,25) linspace(0,1,25)],[fliplr(d13) d13],'LineWidth',1,'Color',[98 12 232]/255)
            plot([linspace(-1,0,25) linspace(0,1,25)],[fliplr(d12) d12],'LineWidth',1,'Color',[12 111 232]/255)
            plot([linspace(-1,0,25) linspace(0,1,25)],[fliplr(d11) d11],'LineWidth',1,'Color',[0 255 255]/255)
            
            plot([linspace(-1,0,51) linspace(0,1,51)],[flipud(n14); n14],'LineWidth',2,'Color',[209 0 255]/255)
            plot([linspace(-1,0,36) linspace(0,1,36)],[flipud(n13); n13],'LineWidth',2,'Color',[98 12 232]/255)
            plot([linspace(-1,0,26) linspace(0,1,26)],[flipud(n12); n12],'LineWidth',2,'Color',[12 111 232]/255)
            plot([linspace(-1,0,19) linspace(0,1,19)],[flipud(n11); n11],'LineWidth',2,'Color',[0 255 255]/255)
            title('Gradient Comparison (Normalized)','FontSize',12,'FontName','Arial')
            xlabel('DV Coordinate ','FontSize',12,'FontName','Arial')
            ylabel('Intensity (AU)','FontSize',12,'FontName','Arial')
            set(gcf,'Name','Gradient_comparison_normalized')
            set(gca,'FontSize',12,'FontName','Arial')
            box on
            
            
            
            s14=protein.dlNuc.NC14(:,end);
            s13=protein.dlNuc.NC13(:,end);
            s12=protein.dlNuc.NC12(:,end);
            s11=protein.dlNuc.NC11(:,end);
            
            figure
            plot([linspace(-1,0,51) linspace(0,1,51)],Beta*[flipud(s14); s14],'LineWidth',2,'Color',[209 0 255]/255)
            hold on
            plot([linspace(-1,0,36) linspace(0,1,36)],Beta*[flipud(s13);s13],'LineWidth',2,'Color',[98 12 232]/255)
            plot([linspace(-1,0,26) linspace(0,1,26)],Beta*[flipud(s12);s12],'LineWidth',2,'Color',[12 111 232]/255)
            plot([linspace(-1,0,19) linspace(0,1,19)],Beta*[flipud(s11);s11],'LineWidth',2,'Color',[0 255 255]/255)
            
            
            d14=protein.dlNuc.NC14(:,end)+protein.dlCactNuc.NC14(:,end);
            d13=protein.dlNuc.NC13(:,end)+protein.dlCactNuc.NC13(:,end);
            d12=protein.dlNuc.NC12(:,end)+protein.dlCactNuc.NC12(:,end);
            d11=protein.dlNuc.NC11(:,end)+protein.dlCactNuc.NC11(:,end);
            
            plot([linspace(-1,0,51) linspace(0,1,51)],Beta*[flipud(d14);d14],'LineWidth',2,'Color',[209 0 255]/255)
            plot([linspace(-1,0,36) linspace(0,1,36)],Beta*[flipud(d13);d13],'LineWidth',2,'Color',[98 12 232]/255)
            plot([linspace(-1,0,26) linspace(0,1,26)],Beta*[flipud(d12);d12],'LineWidth',2,'Color',[12 111 232]/255)
            plot([linspace(-1,0,19) linspace(0,1,19)],Beta*[flipud(d11);d11],'LineWidth',2,'Color',[0 255 255]/255)
            title('Deconvolution of Gradient','FontSize',12,'FontName','Arial')
            xlabel('DV Coordinate ','FontSize',12,'FontName','Arial')
            ylabel('Intensity (AU)','FontSize',12,'FontName','Arial')
            ylim([0 6000])
            set(gcf,'Name','Gradient_deconvolution_in_space')
            set(gca,'FontSize',12,'FontName','Arial')
            box on
            
            % deconvolution of gradient (raw sim)
            figure
            plot([linspace(-1,0,51) linspace(0,1,51)],[flipud(s14); s14],'LineWidth',2,'Color',[209 0 255]/255)
            hold on
            plot([linspace(-1,0,36) linspace(0,1,36)],[flipud(s13);s13],'LineWidth',2,'Color',[98 12 232]/255)
            plot([linspace(-1,0,26) linspace(0,1,26)],[flipud(s12);s12],'LineWidth',2,'Color',[12 111 232]/255)
            plot([linspace(-1,0,19) linspace(0,1,19)],[flipud(s11);s11],'LineWidth',2,'Color',[0 255 255]/255)
            
            
            d14=protein.dlNuc.NC14(:,end)+protein.dlCactNuc.NC14(:,end);
            d13=protein.dlNuc.NC13(:,end)+protein.dlCactNuc.NC13(:,end);
            d12=protein.dlNuc.NC12(:,end)+protein.dlCactNuc.NC12(:,end);
            d11=protein.dlNuc.NC11(:,end)+protein.dlCactNuc.NC11(:,end);
            
            plot([linspace(-1,0,51) linspace(0,1,51)],[flipud(d14);d14],'LineWidth',2,'Color',[209 0 255]/255)
            plot([linspace(-1,0,36) linspace(0,1,36)],[flipud(d13);d13],'LineWidth',2,'Color',[98 12 232]/255)
            plot([linspace(-1,0,26) linspace(0,1,26)],[flipud(d12);d12],'LineWidth',2,'Color',[12 111 232]/255)
            plot([linspace(-1,0,19) linspace(0,1,19)],[flipud(d11);d11],'LineWidth',2,'Color',[0 255 255]/255)
            
            title('Deconvolution of Gradient','FontSize',12,'FontName','Arial')
            xlabel('DV Coordinate ','FontSize',12,'FontName','Arial')
            ylabel('Intensity (AU)','FontSize',12,'FontName','Arial')
            
            set(gcf,'Name','Gradient_deconvolution_in_space (raw sim)')
            set(gca,'FontSize',12,'FontName','Arial')
            box on
            
            % deconvolution of gradient (raw sim) nc14 only
            figure
            plot([linspace(-1,0,51) linspace(0,1,51)],[flipud(s14); s14],'LineWidth',2,'Color',[209 0 255]/255)
            hold on
            
            d14=protein.dlNuc.NC14(:,end)+protein.dlCactNuc.NC14(:,end);
            dd14 = protein.dlCactNuc.NC14(:,end);
            plot([linspace(-1,0,51) linspace(0,1,51)],[flipud(d14);d14],'LineWidth',2,'Color',[209 0 255]/255)
            plot([linspace(-1,0,51) linspace(0,1,51)],[flipud(dd14);dd14],'LineWidth',2,'Color',[209 0 255]/255)
            
            
            title('Deconvolution of Gradient','FontSize',12,'FontName','Arial')
            xlabel('DV Coordinate ','FontSize',12,'FontName','Arial')
            ylabel('Intensity (AU)','FontSize',12,'FontName','Arial')
            
            set(gcf,'Name','Gradient_deconvolution_in_space (raw sim)')
            set(gca,'FontSize',12,'FontName','Arial')
            box on
            
            % if we're talking about gene expression, this will plot the
            % gene activation parameters against the dl gradient
            if ~isempty(params)
                figure
                plot(linspace(0,1,51),s14,'LineWidth',2,'Color',[209 0 255]/255)
                hold on
                %                 plot(linspace(0,1,36),s13,'LineWidth',2,'Color',[98 12 232]/255)
                %                 plot(linspace(0,1,26),s12,'LineWidth',2,'Color',[12 111 232]/255)
                %                 plot(linspace(0,1,19),s11,'LineWidth',2,'Color',[0 255 255]/255)
                
                plot([0,1],[params(1) params(1)],'Color',[153 0 0]/255)
                plot([0,1],[params(2) params(2)],'Color',[0 0 204]/255)
                plot([0,1],[params(3) params(3)],'Color',[0 153 0]/255)
                plot([0,1],[params(4) params(4)],'Color',[204 204 0]/255)
                set(gca,'YScale','log')
                title('Gene expression params','FontSize',12,'FontName','Arial')
                xlabel('DV Coordinate ','FontSize',12,'FontName','Arial')
                ylabel('Intensity (AU)','FontSize',12,'FontName','Arial')
                
                set(gcf,'Name','Gene_exp_vs_gradient_decon')
                set(gca,'ylim',[min(s14) max(s14)])
                
                figure
                plot(linspace(0,1,51),d14,'LineWidth',2,'Color',[209 0 255]/255)
                hold on
                %                 plot(linspace(0,1,36),d13,'LineWidth',2,'Color',[98 12 232]/255)
                %                 plot(linspace(0,1,26),d12,'LineWidth',2,'Color',[12 111 232]/255)
                %                 plot(linspace(0,1,19),d11,'LineWidth',2,'Color',[0 255 255]/255)
                
                plot([0,1],[params(1) params(1)],'Color',[153 0 0]/255)
                plot([0,1],[params(2) params(2)],'Color',[0 0 204]/255)
                plot([0,1],[params(3) params(3)],'Color',[0 153 0]/255)
                plot([0,1],[params(4) params(4)],'Color',[204 204 0]/255)
                set(gca,'YScale','log')
                title('Gene expression params','FontSize',12,'FontName','Arial')
                xlabel('DV Coordinate ','FontSize',12,'FontName','Arial')
                ylabel('Intensity (AU)','FontSize',12,'FontName','Arial')
                
                set(gcf,'Name','Gene_exp_vs_gradient_naive')
                %                 set(gca,'ylim',[min(d14,params) max(d14)])
            end
            
            %
            figure
            relativeError14 = abs((d14-s14)./d14);
            relativeError13 = abs((d13-s13)./d13);
            relativeError12 = abs((d12-s12)./d12);
            relativeError11 = abs((d11-s11)./d11);
            hold on
            plot([linspace(-1,0,51) linspace(0,1,51)],[flipud(relativeError14);relativeError14],'LineWidth',2,'Color',[209 0 255]/255)
            plot([linspace(-1,0,36) linspace(0,1,36)],[flipud(relativeError13);relativeError13],'LineWidth',2,'Color',[98 12 232]/255)
            plot([linspace(-1,0,26) linspace(0,1,26)],[flipud(relativeError12);relativeError12],'LineWidth',2,'Color',[12 111 232]/255)
            plot([linspace(-1,0,19) linspace(0,1,19)],[flipud(relativeError11);relativeError11],'LineWidth',2,'Color',[0 255 255]/255)
            %           plot([linspace(-1,0,51) linspace(0,1,51)],[flipud(relativeError);relativeError])
            title('Relative error in gradient','FontSize',12,'FontName','Arial')
            xlabel('DV Coordinate ','FontSize',12,'FontName','Arial')
            ylabel('Percent error','FontSize',12,'FontName','Arial')
            legend('NC14','NC13','NC12','NC11')
            set(gca,'FontSize',12,'FontName','Arial')
            set(gcf,'Name','Relative_error_in_gradient')
            box on
    end

    function plotGeneExpression3D
        
    end

    function plotProtein3D
        for j = 1:length(ncs)
            if j ==1
                hold off
            end
            surf(time.(ts{j})',time.(xs{j})',Beta*protein.(plot_type).(ncs{j}));
            if j == 1
                hold on
            end
            surf(time.(ts{j})',-time.(xs{j})',Beta*(protein.(plot_type).(ncs{j})));
        end
        shading flat
    end

    function plotProtein2D
        for j = 1:length(ncs)
            plot(time.(ts{j})',protein.(plot_type).(ncs{j}));
        end

    end
end