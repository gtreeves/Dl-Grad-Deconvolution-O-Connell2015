% script_canonical_gene_plot
%
% This script assumes the delta-factor is 1 for each gene.  This may not be
% true if the canonical gene expression profile does not reflect the true
% width of the gene expression pattern.

figure
COLOR = {'r','b','c','g','y'};
genename = {'sna','vnd','ind','sog','dpp'};
s0 = [0 0.27 0.42 0.3 1];
for i = 1:length(COLOR)
	matname = [genename{i},'avg.mat'];
	load(matname)
	plot(s-s_offset+s0(i),t,'Color',COLOR{i},'Linewidth',2)
	hold on
end
xlim([0 1])
ylim([0 1])