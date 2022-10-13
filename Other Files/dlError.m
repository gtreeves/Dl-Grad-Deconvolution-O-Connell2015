function [error,Beta] = dlError(protein)

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
    
    error = sqrt((e11n+e12n+e13n+e14n)/length(D1));
    
catch
    disp('Vectors not the same size')
    error = 1e12;
    Beta = 0;
end

end