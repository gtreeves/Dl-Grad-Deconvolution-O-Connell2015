
load thiele_transitionwidth_stuff
eta = 0.2;
Dn = 0.02;
thiele = abs(dcdx)./c;
figure
plot(x,[c,abs(dcdx)])
hold on
plot(x,(eta./sqrt(c))./(thiele*Dn),'r')
ylim([0 20])
legend('c (NC14)','abs(dcdx)','trans width')

% 
% %-- 5/22/2014 2:39 PM --%
% uiopen('C:\Users\gtreeves\Desktop\Thiele modulus.fig',1)
% X = get(gco,'Xdata')';
% Y = get(gco,'Ydata')';
% dcdx = Y;
% c = get(gco,'Ydata')';
% thiele = abs(dcdx)/sqrt(c);
% hold on
% plot(X,N,'r')
% N = -dcdx./sqrt(c);
% plot(X,N,'r')
% plot(X,1./N,'r')
% figure
% plot(X,1./N,'r')
% N(1)
% dcdx(0)
% dcdx(1)
% ylim([0 5])
% ylim([0 1])
% plot(X,1./N/0.02*0.1,'r')
% ylim([0 10])
% ylim([0 5])
% hold on
% plot(X,1./N/0.02*0.3,'g')
% ylim([0 10])
% ylim([0 12])
% %-- 5/23/2014 2:29 PM --%
% uiopen('C:\Users\gtreeves\Desktop\Thiele modulus.fig',1)
% X = get(gco,'Xdata')';
% Y = get(gco,'Ydata')';
% dcdx = Y;
% c = get(gco,'Ydata')';
% N = -dcdx/sqrt(c);
% hold on
% plot(X,N,'r')
% N = -dcdx./sqrt(c);
% plot(X,N,'r')
% plot(X,1./N,'r')
% figure
% plot(X,1./N,'r')
% N(1)
% dcdx(0)
% dcdx(1)
% ylim([0 5])
% ylim([0 1])
% plot(X,1./N/0.02*0.1,'r')
% ylim([0 10])
% ylim([0 5])
% hold on
% plot(X,1./N/0.02*0.3,'g')
% ylim([0 10])
% ylim([0 12])
% X = get(gco,'Xdata')';
% Y = get(gco,'Ydata')';
% dcdx = Y;
% c = get(gco,'Ydata')';
% N = -dcdx./sqrt(c);
% figure
% plot(X,1./N/0.02*0.1,'r')
% ylim([0 12])
% hold on
% plot(X,1./N/0.02*0.3,'g')
% uiopen('C:\Users\gtreeves\Desktop\x error comparison.fig',1)
% ylim([0 12])
% ylim([0 20])
% uiopen('C:\Users\gtreeves\Desktop\x error comparison.fig',1)
% phi = get(gco,'Ydata')';
% plot(X,phi*0.02,'r')
% figure
% plot(X,phi*0.02,'r')
% xlim([0 0.5])
% xlim([0 `])
% xlim([0 1])
% ylim([0 0.])
% ylim([0 0.5])
% ylim([0 0.8])
% ylim([0 0.2])
% ylim([0 0.3])
% figure
% delta_c = -sqrt(c)./dcdx;
% plot(X,delta_c./c)
% delta_c = sqrt(c);
% plot(X,delta_c./c)
% plot(X,0.1*delta_c./c)
% plot(X,phi*0.02,'r')
% ylim([0 0.3])
% hold on
% plot(X,0.1*delta_c./c)
% ylim([0 0.8])
% Deltax = 0.02;
% figure
% plot(X,phi*Deltax,'r')
% hold on
% plot(X,0.1./sqrt(c))
% ylim([0 0.8])
% figure
% axes
% ylim([0 10])
% ylim([0 5])
% 4/3*pi
% (4/3*pi)^(2/3)
% (4/3*pi)^(1/3)