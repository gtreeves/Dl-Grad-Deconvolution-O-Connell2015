% script_checkdeconresults
%
% This script is to check the deconvolution results that Mike O'Connell
% made. There is one nagging feeling about his results. It seems like he
% didn't carry over the final conditions of mitosis over to the next
% interphase. In essence, interphases 11-14 pick up exactly where the
% previous interphase left off. The calamity which is mitosis doesn't even
% happen to his code.
%
% So, this is double-checking that.


addpath('C:\Users\gtreeves\Documents\Dropbox\Matlab\Dorsal\MDO\deconvolution\Other Files')
load('deconvolutionResults.mat')
i = 1;
[error,penalty,protein,XT] = decon_log_nest_JPattern(xb(i,:));

plotComparison(protein.dlNuc,protein.dlCactNuc,XT);
