function plotAdjusted(dlNuc,dlCactNuc,Alpha,Beta)
% 
% Plotting adjusted data with the simulation
% 
figure
subplot(1,2,1)
surf(time.T11',time.X11',Beta*(dlNuc.DN11+dlCactNuc.DCN11)+Alpha)
hold on
surf(time.T12',time.X12',Beta*(dlNuc.DN12+dlCactNuc.DCN12)+Alpha)
surf(time.T13',time.X13',Beta*(dlNuc.DN13+dlCactNuc.DCN13)+Alpha)
surf(time.T14',time.X14',Beta*(dlNuc.DN14+dlCactNuc.DCN14)+Alpha)
ylabel('space')
xlabel('time')
title('dlNuc+dlCactNuc')
shading flat

subplot(1,2,2)
surf(T11D',X11D',N11D)
hold on
surf(T12D',X12D',N12D)
surf(T13D',X13D',N13D)
surf(T14D',X14D',N14D)
ylabel('space')
xlabel('time')
title('Data')
shading flat