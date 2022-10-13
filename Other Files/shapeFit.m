function [A,B,M,sig] = shapeFit(protein,time)

%% Unpack "protein"

dlNuc = protein{1};
dlCactNuc = protein{3};

t{1} = dlNuc.NC10+dlCactNuc.NC10;
t{2} = dlNuc.NC11+dlCactNuc.NC11;
t{3} = dlNuc.NC12+dlCactNuc.NC12;
t{4} = dlNuc.NC13+dlCactNuc.NC13;
t{5} = dlNuc.NC14+dlCactNuc.NC14;

s{1} = time.X10(1,:);
s{2} = time.X11(1,:);
s{3} = time.X12(1,:);
s{4} = time.X13(1,:);
s{5} = time.X14(1,:);

% preallocate
A = cell(5,1); A{1} = zeros(1,length(t{1})); A{2} = zeros(1,length(t{2}));
A{3} = zeros(1,length(t{3})); A{4} = zeros(1,length(t{4})); 
A{5} = zeros(1,length(t{5}));

B = A; M = A; sig = A; 

for i = 1:5
    [~,b] = size(t{i});
    for j = 1:b
        [A{i}(1,j),B{i}(1,j),M{i}(1,j),sig{i}(1,j)] = fitgauss(s{i}',t{i}(:,j));
        
    end
end