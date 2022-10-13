function [dlerror,geneerror,Beta] = totalError(protein)

try
    % Greg's data
    [C1,dC1] = gregsData('nuclear');
    
    % Combine simulation data
    N11S = protein.dlNuc.NC11 + protein.dlCactNuc.NC11;
    N12S = protein.dlNuc.NC12 + protein.dlCactNuc.NC12;
    N13S = protein.dlNuc.NC13 + protein.dlCactNuc.NC13;
    N14S = protein.dlNuc.NC14 + protein.dlCactNuc.NC14;
    
    % individual NCs as columns
    N11SNcol=N11S(:);
    N12SNcol=N12S(:);
    N13SNcol=N13S(:);
    N14SNcol=N14S(:);
    
    % Simulation data as a column vector
    D1 = [N11SNcol;N12SNcol;N13SNcol;N14SNcol];
    
    % nuclear data by nc
    C11Ncol=C1(1:304);
    C12Ncol=C1(305:980);
    C13Ncol=C1(981:2780);
    C14Ncol=C1(2781:end);
    dC11Ncol=dC1(1:304);
    dC12Ncol=dC1(305:980);
    dC13Ncol=dC1(981:2780);
    dC14Ncol=dC1(2781:end);
    
    Beta = mean(D1.*C1)/mean(D1.^2);
    
    e11n = sum(((C11Ncol - Beta*N11SNcol)./dC11Ncol).^2);
    e12n = sum(((C12Ncol - Beta*N12SNcol)./dC12Ncol).^2);
    e13n = sum(((C13Ncol - Beta*N13SNcol)./dC13Ncol).^2);
    e14n = sum(((C14Ncol - Beta*N14SNcol)./dC14Ncol).^2);
    
    dlerror = sqrt((e11n+e12n+e13n+e14n)/length(D1));
    
catch
    disp('dl/Cact error not calculated')
    dlerror = 1e12;
    Beta = 0;
end

try
    
    sim = [protein.sna.NC14(:,end) protein.sog.NC14(:,end) protein.dpp.NC14(:,end) protein.vnd.NC14(:,end)];
    
    sim(isnan(sim))= 0;
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
    
    
    
    
    
    
    % interpolating down to the size of the simulation data
    x1 = linspace(0,1,151);
    x2 = linspace(0,1,51);
    geneavg = interp1(x1,geneavg,x2);
    geneerror = sum((geneavg(:)-sim(:)).^2);
catch
    disp('Gene error not calculated')
    geneerror = 1e12;
end

end