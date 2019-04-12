close all; clear all; clc

dt = 0.6562;   % [in] throat diameter
rt = dt/2;  % [in] throat radius

di = 1.68;  % [in] nozzle inlet diameter
ri = di/2;  % [in] nozzle inlet radius

% converging section
Ld = 1;     % L/D where L is length of converging section and D is inlet diameter
L = Ld*di;
Area_ratio = (ri^2)/(rt^2);
x1 = [0:0.025:L+0.05];
a0 = sqrt(Area_ratio) * rt - rt; 
a1 = 0;
a2 = 0;
a3 = -639875/592704; %-321875/296352;
a4 = 15996875/16595712;
a5 = -79984375/348509952;

a = [a5, a4, a3, a2, a1, a0];
y1 = polyval(a, x1);%a0 + a3*x.^3 + a4*x.^4 + a5*x.^5 + rt;

%plot
plot(x1,y1 + rt,'b'); xlabel('position [in]'); ylabel('radius from centerline [in]');
hold on;
plot(x1,-(y1 + rt),'b');
plotfixer
 
% L = 1.68, 0.5*(1.68 - 0.6562)+ a*L^3 + b*L^4 + c*L^5 = 0; 6*a + 12*b*L+ 20*c*L^2 = 0, 3*a + 4*b*L + 5*c*L^2 = 0

% L = 1.68, 
% 0.3*((1.68/0.6)^2 - 1) + a*L^3 + b*L^4 + c*L^5 = 0; 
% 6*a + 12*b*L+ 20*c*L^2 = 0, 
% 3*a + 4*b*L + 5*c*L^2 = 0


%% diverging section
epsilon = 2.2;  % expansion ratio
At = pi*rt^2;   % area at throat
Ae = epsilon * At;  % area at exit
re = sqrt(Ae/pi);

% L = 1.32; (sqrt(2.25) - 1)*(0.3) = L^3*(a + b*L + c*L^2); 0 = 3a + 6b*L + 10c*L^2; 0 = 3a+4b*L+5c*L^2

L2 = 1.32;          % [in]
% b0 = 0; b1 = 0; b2 = 0; 
% b3 = 0.652183; 
% b4 = -0.741117;
% b5 = 0.224581; 
% 
% b = [b5, b4, b3, b2, b1, b0];
% x2 = linspace(0, L2);
% 
% y2 = polyval(b, x2);

% plot(x2+L, y2 + rt, 'b')
% plot(x2+L, -(y2 + rt), 'b')

plot([1.68 ,3], [rt, re], 'b')
plot([1.68 ,3], -[rt, re], 'b')

plotfixer();
set(gcf,'Units','inches','position',[0 0 9 6])

nozzle = [x1, [1.68 ,3];  y1+rt, [rt, re]];
coordinates = y1*2;

