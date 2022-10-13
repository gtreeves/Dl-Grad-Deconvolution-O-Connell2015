% script_signaltonoise

eta0 = 0.1;
eta1 = 0.2;
ETA0 = [0 0.01 0.03 0.1 0.3];

close all
load signal_to_noise % contains c, dcdx, Dn, L, and x.
v = M < 0;
S = -1./c.*dcdx.*repmat(Dn,1,length(x))./repmat(L,1,length(x)); % signal

% BS
c_bs = c - repmat(B+M,1,length(x));
S_bs = -1./c_bs.*dcdx.*repmat(Dn,1,length(x))./repmat(L,1,length(x));

%
% Color definition
%
C = [1 0 1];
C1 = [0 1 0];
[~,COLR] = colorinterpolate(length(ETA0));
COLR_bs = COLR + repmat(C1,length(ETA0),1); 
COLR = COLR + repmat(C,length(ETA0),1);

% [~,COLR] = colorinterpolate(length(ETA0),'spring');
% [~,COLR_bs] = colorinterpolate(length(ETA0),'winter');
for i = 1:length(ETA0)
	eta0 = ETA0(i);
	
	%
	% Noise calc without bs
	%
	eta = eta1*sqrt(trapz(x',c')'); % an approximation to properly scale eta's
	eta_in = eta0*trapz(x',c')';
	N_ex = repmat(eta,1,length(x))./sqrt(c);
	N_in = repmat(eta_in,1,length(x))./c;
	N = N_in + N_ex;
	
	%
	% Noise calc with bs
	%
	eta = eta1*sqrt(trapz(x',c_bs')'); % an approximation to properly scale eta
	eta_in = eta0*trapz(x',c_bs')';
	N_ex = repmat(eta,1,length(x))./sqrt(c_bs); % this is not properly scaled.
	N_in = repmat(eta_in,1,length(x))./c_bs;
	N_bs = N_in + N_ex;
	
	G = N./S; % number of nuclei in your transition width.
	G_totdl = mean(G(v,:));
	
	G1 = N_bs./S_bs; % number of nuclei in your transition width.
	G_bs = mean(G1(v,:));
	
	h = plot(x,[G_totdl;G_bs],'Linewidth',2);
	hold on
	colr = COLR(i,:); colr = colr/max(colr);
	colr_bs = COLR_bs(i,:); colr_bs = colr_bs/max(colr_bs);
	set(h(1),'Color',colr)
	set(h(2),'Color',colr_bs)
	xlim([0 1])
	ylim([0 30])
	set(gca,'Fontsize',24,'YTick',0:5:30,'XTick',0:0.2:1)
	xlabel('DV coordinate')
	ylabel('Number of nuclei')
	legend('Total dl','Bkgrnd Subtr','Location','Northwest')
	
% 	s = sprintf('\\eta_0 = %g',eta0);
% 	title(s)
end
print(gcf,'effect_of_intrinsic_noise.jpg','-djpeg','-r300')

% figure
% plot(x,G(v,:)) % # of nuclei
% hold on
% G_totdl = mean(G(v,:));
% plot(x,G_totdl,'k','Linewidth',2)
% ylim([0 30])
% xlabel('DV coordinate')
% ylabel('Number of nuclei')
% title('Sensitivity of the Dorsal gradient at each DV point')
% 
% figure
% plot(x,1./G(v,:))
% xlabel('DV coordinate')
% ylabel('signal to noise ratio')
% hold on
% plot(x,1./G_totdl,'k','Linewidth',2)

%% Now do the same thing but with background subtraction as an approx




% figure
% plot(x,G(v,:)) % # of nuclei
% hold on
% G_bs = mean(G(v,:));
% plot(x,G_bs,'k','Linewidth',2)
% ylim([0 30])
% xlabel('DV coordinate')
% ylabel('Number of nuclei - Bg Sub')
% title('Sensitivity of the Dorsal gradient at each DV point')
% 
% figure
% plot(x,1./G(v,:))
% ylim([0 2])
% xlabel('DV coordinate')
% ylabel('signal to noise ratio - Bg Sub')
% hold on
% plot(x,1./G_bs,'k','Linewidth',2)

%% Do it all over again but plot only the means

% figure
% plot(x,[G_totdl;G_bs],'Linewidth',2)
% ylim([0 30])
% xlabel('DV coordinate')
% ylabel('Number of nuclei')
% legend('Total dl','Bkgrnd Subtr')
% title('Sensitivity of the Dorsal gradient at each DV point')
% 
% figure
% plot(x,1./[G_totdl;G_bs],'Linewidth',2)
% ylim([0 2])
% xlabel('DV coordinate')
% ylabel('signal to noise ratio - Bg Sub')
% legend('Total dl','Bkgrnd Subtr')





