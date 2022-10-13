function save_all_figures_to_publication(appendix)
pth = cd;
mfolder = strfind(pth,'MATLAB')-1;
dir_name = strcat(pth(1:mfolder),'Documents/Papers/Publications');
figlist=findobj('type','figure');


for i=1:numel(figlist)
    if strcmp(get(figlist(i),'WindowStyle'),'docked')
        set(figlist(i),'WindowStyle','normal')
    end
    ratio=1;
    set(figlist(i),'Position',[0 50 ratio*250 ratio*195])
    set(figlist(i),'Paperpositionmode','auto')
    figname = get(figlist(i),'Name');
    if isempty(figname)
%     saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '_' appendix '.fig']));
    saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '_' appendix '.eps']),'epsc')
    else
%     saveas(figlist(i),fullfile(dir_name,['figure_' figname '_' appendix '.fig']));
    saveas(figlist(i),fullfile(dir_name,['figure_' figname '_' appendix '.eps']),'epsc')        
    end
end
end