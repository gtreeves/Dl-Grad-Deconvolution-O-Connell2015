function [error,varargout] = mRNAerrorF(protein)

t = [];

n = 4;

genename = {'sna', 'sog', 'dpp', 'vnd'};
% s0 = [0 0.3 1 0.27]; % from Greg's code

% initialize values for loop
geneavg = zeros(151,n);

for i = 1:n
    matname = [genename{i},'avg.mat'];
    load(matname)
    if i == 2
        geneavg(1:121,i) = t(31:end);
    else
        geneavg(:,i) = t;
    end
    
end

geneavg(1:25,1)=round(geneavg(1:25,1));
geneavg(133:151,3)=round(geneavg(133:151,3));



sim = [protein.sna.NC14(:,end) protein.sog.NC14(:,end)...
    protein.dpp.NC14(:,end) protein.vnd.NC14(:,end)];

sim(isnan(sim))= 0;


% interpolating down to the size of the simulation data
x1 = linspace(0,1,151);
x2 = linspace(0,1,51);
geneavg = interp1(x1,geneavg,x2);



error = sum((geneavg(:)-sim(:)).^2);

if nargout>1
    varargout = cell(1,nargout-1);
    argsOut = {x2 geneavg x1 sim};
    for i = 1:nargout-1
        varargout{i}=argsOut{i};
    end
end
end