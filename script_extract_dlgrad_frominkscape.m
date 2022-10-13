% script_extract_dlgrad_frominkscape
%
% This script takes the coordinates of the free dl gradient from an
% inkscape image (see Fig 5 from O'Connell and Reeves, 2015) and converts
% them to plot coordinates in matlab.

%
% Initial coordinates
%
p = [35.895448,348.12205 2.201953,0.21424 2.202295,0.42846 2.201954,0.71383 2.201953,1.00005 2.202295,1.2854 2.201954,1.57077 2.201953,1.78585 2.202295,2.07036 2.201954,2.28545 2.201953,2.49969 2.202295,2.64193 2.201953,2.78505 2.201954,2.92815 2.202295,2.99911 2.201953,3.07057 2.201954,3.14222 2.202295,3.14187 2.201953,3.07092 2.201954,3.14187 2.202295,3.07092 2.201953,2.99928 2.201954,2.92764 2.202294,2.85668 2.201954,2.78504 2.273386,2.64211 2.201783,2.57081 2.201783,2.49951 2.202636,2.28494 2.201783,2.14234 2.201786,2.07105 2.20264,1.85681 2.20178,1.78517 2.20178,1.57093 2.20264,1.42835 2.20178,1.2854 2.20178,1.07117 2.20264,0.99988 2.20178,0.85693 2.20179,0.7853 2.20263,0.6427 2.20178,0.50011 2.20179,0.42847 2.20263,0.42847 2.20179,0.28553 2.20178,0.21423 2.20264,0.21424 2.20178,0.14293 2.20178,0.0713 2.20264,0.0713 2.27253,0.0716];

%
% Extract x and y values
%
x = p(1:2:end);
y = p(2:2:end);

%
% It turns out the coordinates here are formatted the following way. The
% first coordinate is location, and the remaining ones are displacements
% from the previous point. So need a cumulative sum.
%
x = cumsum(x);
y = cumsum(y);

%
% Now, coordinates for the plot itself.
%
t1 = [35.895448,414.74888 3.052772,0]; % tickmark at x = 0, y = 10^-1 
t2 = [146.1406,374.75839 -3.375,0]; % tickmark at x = 1, y = 10^0

x0 = t1(1); ym1 = t1(2);
x1 = t2(1); y0 = t2(2);

y = 10.^(1/(y0 - ym1)*(y - y0)); y = y/max(y);
x = 1/(x1 - x0)*(x - x0);

save dlgrad_deconv_Fig5 x y

% figure
% plot(xnew,ynew)
% xlim([0 1])
% ylim([3e-2 5e0])







