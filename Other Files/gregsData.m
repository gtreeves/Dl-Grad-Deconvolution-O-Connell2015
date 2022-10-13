function [C,dC] = gregsData(options,varargin)

% gregsData(options) - a function that creates a column vector of Greg's
% data useful for calculating error
%
% C : the column vector of data
% dC: the weights of each data point
%
% options: 'nuclear' - nuclear dl fluorescence
%          'totalDl' - total dl fluorescence

load('hybrid_embryo.mat')
switch options
    case 'nuclear'
        A = data.A;
        B = data.B;
        dA = data.dA;
        dB = data.dB;
    case {'totalDl','totalDorsal','total'}
        A = data.A2;
        B = data.B2;
        dA = data.dA;
        dB = data.dB;
end

args = struct('M',[19 26 36 51],'yesplot',false);
if nargin > 1
    if mod(length(varargin),2)==0
        for i = 1:2:length(varargin)
            args.(varargin{i})=varargin{i+1};
        end
    else
        error('Varargin requires name-value pairs. Some options not specified')
    end
end
%% NC 11
% N11T = t(1:16);
N11A = A(1:16);
N11B = B(1:16);
dA11 = dA(1:16);
dB11 = dB(1:16);

sig = 0.15;
m = -0.07;
M = args.M(1);
x = linspace(0,1,M);
N11D = repmat(N11A,1,M).*exp(-repmat(x,length(N11A),1).^2/2/sig^2) + ...
    repmat(N11B,1,M) + m*N11A*x;
N11D = fliplr(rot90(N11D,3));

dC11 = exp(-repmat(x,length(N11A),1).^2/2/sig^2).*repmat(dA11,1,length(x))+m*dA11*x + m*repmat(x,length(dA11),1) + repmat(dB11,1,length(x));

%% NC 12
% N12T = t(34:59);
N12A = A(34:59);
N12B = B(34:59);
dA12 = dA(34:59);
dB12 = dB(34:59);
M = args.M(2);
x = linspace(0,1,M);
N12D = repmat(N12A,1,M).*exp(-repmat(x,length(N12A),1).^2/2/sig^2) + ...
    repmat(N12B,1,M) + m*N12A*x;
N12D = fliplr(rot90(N12D,3));
dC12 = exp(-repmat(x,length(N12A),1).^2/2/sig^2).*repmat(dA12,1,length(x))+m*dA12*x + m*repmat(x,length(dA12),1) + repmat(dB12,1,length(x));

%% NC 13
% N13T = t(77:126);
N13A = A(77:126);
N13B = B(77:126);
dA13 = dA(77:126);
dB13 = dB(77:126);
M = args.M(3);
x = linspace(0,1,M);
N13D = repmat(N13A,1,M).*exp(-repmat(x,length(N13A),1).^2/2/sig^2) + ...
    repmat(N13B,1,M) + m*N13A*x;
N13D = fliplr(rot90(N13D,3));
dC13 = exp(-repmat(x,length(N13A),1).^2/2/sig^2).*repmat(dA13,1,length(x))+m*dA13*x + m*repmat(x,length(dA13),1) + repmat(dB13,1,length(x));

%% NC 14
% N14T = t(149:332);
N14A = A(149:332);
N14B = B(149:332);
dA14 = dA(149:332);
dB14 = dB(149:332);

M = args.M(4);
x = linspace(0,1,M);
N14D = repmat(N14A,1,M).*exp(-repmat(x,length(N14A),1).^2/2/sig^2) + ...
    repmat(N14B,1,M) + m*N14A*x;
N14D = fliplr(rot90(N14D,3));
dC14 = exp(-repmat(x,length(N14A),1).^2/2/sig^2).*repmat(dA14,1,length(x))+m*dA14*x + m*repmat(x,length(dA14),1) + repmat(dB14,1,length(x));

% Fluorescence data as a column vector
C = [N11D(:);N12D(:);N13D(:);N14D(:)];
dC = [dC11(:);dC12(:);dC13(:);dC14(:)];
% totalDorsal_int = [];

if args.yesplot
    
end
end







% switch options
%     case 'nuclear'
% A = data.A;
% B = data.B;
% dA = data.dA;
% dB = data.dB;
% 
% %% NC 11
% % N11T = t(1:16);
% N11A = A(1:16);
% N11B = B(1:16);
% dA11 = dA(1:16);
% dB11 = dB(1:16);
% 
% sig = 0.15;
% m = -0.07;
% M = 19;
% x = linspace(0,1,M);
% N11D = repmat(N11A,1,M).*exp(-repmat(x,length(N11A),1).^2/2/sig^2) + ...
%     repmat(N11B,1,M) + m*N11A*x;
% N11D = fliplr(rot90(N11D,3));
% 
% dC11 = exp(-repmat(x,length(N11A),1).^2/2/sig^2).*repmat(dA11,1,length(x))+m*dA11*x + m*repmat(x,length(dA11),1) + repmat(dB11,1,length(x));
% 
% %% NC 12
% % N12T = t(34:59);
% N12A = A(34:59);
% N12B = B(34:59);
% dA12 = dA(34:59);
% dB12 = dB(34:59);
% M = 26;
% x = linspace(0,1,M);
% N12D = repmat(N12A,1,M).*exp(-repmat(x,length(N12A),1).^2/2/sig^2) + ...
%     repmat(N12B,1,M) + m*N12A*x;
% N12D = fliplr(rot90(N12D,3));
% dC12 = exp(-repmat(x,length(N12A),1).^2/2/sig^2).*repmat(dA12,1,length(x))+m*dA12*x + m*repmat(x,length(dA12),1) + repmat(dB12,1,length(x));
% 
% %% NC 13
% % N13T = t(77:126);
% N13A = A(77:126);
% N13B = B(77:126);
% dA13 = dA(77:126);
% dB13 = dB(77:126);
% M = 36;
% x = linspace(0,1,M);
% N13D = repmat(N13A,1,M).*exp(-repmat(x,length(N13A),1).^2/2/sig^2) + ...
%     repmat(N13B,1,M) + m*N13A*x;
% N13D = fliplr(rot90(N13D,3));
% dC13 = exp(-repmat(x,length(N13A),1).^2/2/sig^2).*repmat(dA13,1,length(x))+m*dA13*x + m*repmat(x,length(dA13),1) + repmat(dB13,1,length(x));
% 
% %% NC 14
% % N14T = t(149:332);
% N14A = A(149:332);
% N14B = B(149:332);
% dA14 = dA(149:332);
% dB14 = dB(149:332);
% 
% M = 51;
% x = linspace(0,1,M);
% N14D = repmat(N14A,1,M).*exp(-repmat(x,length(N14A),1).^2/2/sig^2) + ...
%     repmat(N14B,1,M) + m*N14A*x;
% N14D = fliplr(rot90(N14D,3));
% dC14 = exp(-repmat(x,length(N14A),1).^2/2/sig^2).*repmat(dA14,1,length(x))+m*dA14*x + m*repmat(x,length(dA14),1) + repmat(dB14,1,length(x));
% 
% % Fluorescence data as a column vector
% C = [N11D(:);N12D(:);N13D(:);N14D(:)];
% dC = [dC11(:);dC12(:);dC13(:);dC14(:)];
% totalDorsal_int = [];
% 
%     case 'totalDl'
%      A2 = data.A2;
% B2 = data.B2;
% sig = data.Sig2;
% 
% % Keeping this as ones, so if we want it later, we still have it in there.
% dA = ones(size(A2));
% dB = dA;
% 
% 
% %% NC 11 (interphase & mitosis)
% % N11T = t(1:33);
% N11A = A2(1:33);
% N11B = B2(1:33);
% dA11 = dA(1:33);
% dB11 = dB(1:33);
% sig11 = sig(1:33);
% 
% m = -0.07;
% M = 19;
% x = linspace(0,1,M);
% N11D = repmat(N11A,1,M).*exp(-repmat(x,length(N11A),1).^2/2./repmat(sig11,1,M).^2) + ...
%     repmat(N11B,1,M) + m*N11A*x;
% N11D = fliplr(rot90(N11D,3));
% 
% dC11 = exp(-repmat(x,length(N11A),1).^2/2./repmat(sig11,1,M).^2).*repmat(dA11,1,length(x))...
%     +m*dA11*x + m*repmat(x,length(dA11),1) + repmat(dB11,1,length(x));
% 
% integral11D = trapz(x,N11D);
% 
% 
% %% NC 12 (interphase & mitosis)
% % N12T = t(34:76);
% N12A = A2(34:76);
% N12B = B2(34:76);
% dA12 = dA(34:76);
% dB12 = dB(34:76);
% sig12 = sig(34:76);
% % N12T = t(33:76);
% % N12A = A2(33:76);
% % N12B = B2(33:76);
% % dA12 = dA(33:76);
% % dB12 = dB(33:76);
% % sig12 = sig(33:76);
% 
% M = 26;
% x = linspace(0,1,M);
% N12D = repmat(N12A,1,M).*exp(-repmat(x,length(N12A),1).^2/2./repmat(sig12...
%     ,1,M).^2) + repmat(N12B,1,M) + m*N12A*x;
% N12D = fliplr(rot90(N12D,3));
% dC12 = exp(-repmat(x,length(N12A),1).^2/2./repmat(sig12,1,M).^2).*...
%     repmat(dA12,1,length(x))+m*dA12*x + m*repmat(x,length(dA12),1)...
%     + repmat(dB12,1,length(x));
% 
% integral12D = trapz(x,N12D);
% 
% 
% %% NC 13 (interphase & mitosis)
% % N13T = t(77:148);
% N13A = A2(77:148);
% N13B = B2(77:148);
% dA13 = dA(77:148);
% dB13 = dB(77:148);
% sig13 = sig(77:148);
% % N13T = t(76:148);
% % N13A = A2(76:148);
% % N13B = B2(76:148);
% % dA13 = dA(76:148);
% % dB13 = dB(76:148);
% % sig13 = sig(76:148);
% 
% M = 36;
% x = linspace(0,1,M);
% N13D = repmat(N13A,1,M).*exp(-repmat(x,length(N13A),1).^2/2./repmat(sig13,1,M).^2) + ...
%     repmat(N13B,1,M) + m*N13A*x;
% N13D = fliplr(rot90(N13D,3));
% dC13 = exp(-repmat(x,length(N13A),1).^2/2./repmat(sig13,1,M).^2).*repmat(dA13,1,length(x))+m*dA13*x + m*repmat(x,length(dA13),1) + repmat(dB13,1,length(x));
% 
% integral13D = trapz(x,N13D);
% 
% 
% %% NC 14 (interphase & mitosis)
% % N14T = t(149:332);
% N14A = A2(149:332);
% N14B = B2(149:332);
% dA14 = dA(149:332);
% dB14 = dB(149:332);
% sig14 = sig(149:332);
% 
% % N14T = t(148:end);
% % N14A = A2(148:end);
% % N14B = B2(148:end);
% % dA14 = dA(148:end);
% % dB14 = dB(148:end);
% % sig14 = sig(148:end);
% 
% M = 51;
% x = linspace(0,1,M);
% N14D = repmat(N14A,1,M).*exp(-repmat(x,length(N14A),1).^2/2./repmat(sig14,1,M).^2) + ...
%     repmat(N14B,1,M) + m*N14A*x;
% N14D = fliplr(rot90(N14D,3));
% dC14 = exp(-repmat(x,length(N14A),1).^2/2./repmat(sig14,1,M).^2).*...
%     repmat(dA14,1,length(x))+m*dA14*x + m*repmat(x,length(dA14),1)...
%     + repmat(dB14,1,length(x));
% 
% integral14D = trapz(x,N14D);
% 
% totalDorsal_int = [integral11D integral12D integral13D integral14D];
% 
% % Fluorescence data as a column vector
% C = [N11D(:);N12D(:);N13D(:);N14D(:)];
% dC = [dC11(:);dC12(:);dC13(:);dC14(:)];   
%     case 'smooth'
%         A = data.A;
% B = data.B;
% dA = data.dA;
% dB = data.dB;
% 
% %% NC 11
% % N11T = t(1:16);
% N11A = A(1:16);
% N11B = B(1:16);
% dA11 = dA(1:16);
% dB11 = dB(1:16);
% 
% sig = 0.15;
% m = -0.07;
% M = 19;
% x = linspace(0,1,M);
% N11A_ = polyfit(linspace(0,1,length(N11A))',N11A,2);
% N11A = polyval(N11A_,linspace(0,1,length(N11A))');
% N11D = repmat(N11A,1,M).*exp(-repmat(x,length(N11A),1).^2/2/sig^2) + ...
%     repmat(N11B,1,M) + m*N11A*x;
% N11D = fliplr(rot90(N11D,3));
% 
% dC11 = ones(size(N11D));
% 
% %% NC 12
% % N12T = t(34:59);
% N12A = A(34:59);
% N12B = B(34:59);
% dA12 = dA(34:59);
% dB12 = dB(34:59);
% M = 26;
% x = linspace(0,1,M);
% N12A_ = polyfit(linspace(0,1,length(N12A))',N12A,2);
% N12A = polyval(N12A_,linspace(0,1,length(N12A))');
% N12D = repmat(N12A,1,M).*exp(-repmat(x,length(N12A),1).^2/2/sig^2) + ...
%     repmat(N12B,1,M) + m*N12A*x;
% N12D = fliplr(rot90(N12D,3));
% dC12 = ones(size(N12D));
% 
% %% NC 13
% % N13T = t(77:126);
% N13A = A(77:126);
% N13B = B(77:126);
% dA13 = dA(77:126);
% dB13 = dB(77:126);
% M = 36;
% x = linspace(0,1,M);
% N13A_ = polyfit(linspace(0,1,length(N13A))',N13A,2);
% N13A = polyval(N13A_,linspace(0,1,length(N13A))');
% N13D = repmat(N13A,1,M).*exp(-repmat(x,length(N13A),1).^2/2/sig^2) + ...
%     repmat(N13B,1,M) + m*N13A*x;
% N13D = fliplr(rot90(N13D,3));
% dC13 = ones(size(N13D));
% 
% %% NC 14
% % N14T = t(149:332);
% N14A = A(149:332);
% N14B = B(149:332);
% dA14 = dA(149:332);
% dB14 = dB(149:332);
% 
% M = 51;
% x = linspace(0,1,M);
% N14A_ = polyfit(linspace(0,1,length(N14A))',N14A,2);
% N14A = polyval(N14A_,linspace(0,1,length(N14A))');
% N14D = repmat(N14A,1,M).*exp(-repmat(x,length(N14A),1).^2/2/sig^2) + ...
%     repmat(N14B,1,M) + m*N14A*x;
% N14D = fliplr(rot90(N14D,3));
% dC14 = ones(size(N14D));
% 
% % Fluorescence data as a column vector
% C = [N11D(:);N12D(:);N13D(:);N14D(:)];
% dC = [dC11(:);dC12(:);dC13(:);dC14(:)];
% totalDorsal_int = [];
% end

% end
