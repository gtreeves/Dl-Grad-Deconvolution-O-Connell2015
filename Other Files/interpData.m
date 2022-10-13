function D = interpData(dlNuc)
    %% Error between 3D surfaces:

    %
    % Interpolate to match data
    %
        N11i = zeros(size(N11D));
        for j = 1:length(dlNuc.DN11(:,1))
            x = linspace(0,1,length(N11D));
            N11i(j,:) = interp1(linspace(0,1,length(dlNuc.DN11(j,:))),dlNuc.DN11(j,:)+dlCactNuc.DCN11(j,:),x);
        end

        N12i = zeros(size(N12D));
        for j = 1:length(dlNuc.DN12(:,1))
            x = linspace(0,1,length(N12D));
            N12i(j,:) = interp1(linspace(0,1,length(dlNuc.DN12(j,:))),dlNuc.DN12(j,:)+dlCactNuc.DCN12(j,:),x);
        end

        N13i = zeros(size(N13D));
        for j = 1:length(dlNuc.DN13(:,1))
            x = linspace(0,1,length(N13D));
            N13i(j,:) = interp1(linspace(0,1,length(dlNuc.DN13(j,:))),dlNuc.DN13(j,:)+dlCactNuc.DCN13(j,:),x);
        end

        N14i = zeros(size(N14D));
        for j = 1:length(dlNuc.DN14(:,1))
            x = linspace(0,1,length(N14D));
            N14i(j,:) = interp1(linspace(0,1,length(dlNuc.DN14(j,:))),dlNuc.DN14(j,:)+dlCactNuc.DCN14(j,:),x);
        end
    
        % interpolated simulation data as a column vector
        D = [N11i(:);N12i(:);N13i(:);N14i(:)]; 
        
end