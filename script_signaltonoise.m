% script_signaltonoise

eta1 = 0.2;
chi = [1 0.99 0.95 0.9 0.7 0.5]; % amt of bg to subtract

close all
load signal_to_noise % contains c, dcdx, Dn, L, and x.
v = M < 0;

%
% Signal without bs
%
S = -1./c.*dcdx.*repmat(Dn,1,length(x))./repmat(L,1,length(x)); % signal

%
% Noise calc without bs
%
eta = eta1*sqrt(trapz(x',c')'); % an approximation to properly scale eta's
N = repmat(eta,1,length(x))./sqrt(c);
G = N./S; % number of nuclei in your transition width.
G_totdl = mean(G(v,:));

%
% Color definition
%
C = [1 0 1];
C1 = [0 1 0];
[~,COLR] = colorinterpolate(length(chi));
COLR_bs = COLR + repmat(C1,length(chi),1); 
COLR = COLR + repmat(C,length(chi),1);

% [~,COLR] = colorinterpolate(length(ETA0),'spring');
% [~,COLR_bs] = colorinterpolate(length(ETA0),'winter');

%
% Plot total dl
%
h = plot(x,G_totdl,'Linewidth',2);
hold on
colr = COLR(i,:); colr = colr/max(colr);
set(h(1),'Color',colr)


%
% Background
%
BG = repmat(B+M,1,length(x));

for i = 1:length(chi)	
	
	%
	% Signal calc with bs
	%
	c_bs = c - chi(i)*BG;
	S_bs = -1./c_bs.*dcdx.*repmat(Dn,1,length(x))./repmat(L,1,length(x));
	
	%
	% Noise calc with bs
	%
	eta = eta1*sqrt(trapz(x',c_bs')'); % an approximation to properly scale eta
	N_bs = repmat(eta,1,length(x))./sqrt(c_bs); % this is not properly scaled.
	G1 = N_bs./S_bs; % number of nuclei in your transition width.
	G_bs = mean(G1(v,:));
	
	h = plot(x,G_bs,'Linewidth',2);
	hold on
	colr_bs = COLR_bs(i,:); colr_bs = colr_bs/max(colr_bs);
	set(h,'Color',colr_bs)
	
% 	s = sprintf('\\eta_0 = %g',eta0);
% 	title(s)
end
xlim([0 1])
ylim([0 100])
% set(gca,'Fontsize',24,'YTick',0:5:30,'XTick',0:0.2:1)
xlabel('DV coordinate')
ylabel('Number of nuclei')
legend('Total dl','Bkgrnd Subtr','Location','Northwest')
print(gcf,'effect_of_bgsub.jpg','-djpeg','-r300')





