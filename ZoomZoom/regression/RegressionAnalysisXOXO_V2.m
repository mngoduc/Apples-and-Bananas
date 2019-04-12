% FUEL GRAIN REGRESSION ANALYSIS
% Computes hole size on fuel grain over time.  Only works for circular holes.
%

% Plots:
% 1)Fuel burned vs. time
% 2)Port perimeter vs. time
% 3)Total port area vs. time
% 4)Mixture ratio vs. time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; 

%init constants
C.burn_time = 3.4;    	% PLACEHOLDER NEED TO BE FIXED
C.time_step = 0.1;        % PLACEHOLDER NEED TO BE FIXED
pause_time = 0.025;
inches_to_m = 0.0254;
bw_plot = 1;
% C. for constants 
C.psi2Pa = 6894.7;      % [psi] -> [Pa]
C.Patm = 101325;        % [Pa]
C.g_rad = 1.005;        % grain outer radius, [in]

% Pull data on the grain design from .mat files
% T. for test
T.filename = input('Enter .mat file of design: '); 
T.FG_data = load(T.filename,'-mat');

C.mult = linspace(0.75,1.5,50);
C.in2m = 0.0254;
C.length = 5.375 * C.in2m;

for j = 1:length(C.mult)
T.FG_data.ri = T.FG_data.r * C.mult(j);
T.hole_data(:,1) = T.FG_data.x;
T.hole_data(:,2) = T.FG_data.y;
T.hole_data(:,3) = T.FG_data.ri;
T.num_holes = length(T.hole_data);

% Start combustion
% D. for test fire data -- incorporate test data 
D = load('testfire7.mat');

T.time = 0:C.time_step:C.burn_time;

T.chamP_abs = (C.Patm+(D.chamP.*C.psi2Pa)).*1E-6; %MPa - lab

% get average chamber pressure over firing

D.thrust(D.thrust < 5) = 0;
for i=1:length(D.thrust)
    if D.thrust(i) == 0
        T.chamP_abs(i) = 0;
    end
end
T.chamP_abs(T.chamP_abs == 0) = [];
% Chamber pressure assumed to be mean value from lab data, can be changed
% for higher accuracy.
T.P = ones(size(T.time)).*mean(T.chamP_abs);  %MPa  

% mdotO2 assumed to be mean value from lab data. 
T.mdotO2 = ones(size(T.time))*mean(D.m_dot_O2);  %kg/s

figure(bw_plot);
hold on;

% Draw hole at each time interval with pause time of 0.5 sec
% use equation r=a(G)^n [cm/s] where G=(mdot_O2/Aport) [kg/m^2s]
% coefficients interpolated over lab data for higher accuracy
filename = 'regressionRate.mat';%input('Enter .mat file of test: '); 
F = load(filename,'-mat');

%convert inches to cm
C.p_HDPE = 950; %kg/m^4
C.FG_length = 5.350 * (2.54e-2); %m

%init arrays
T.r(:,1) = T.FG_data.ri;  
T.m_burned = zeros(size(T.time));
T.Aport = zeros(size(T.time));
T.Aport_inches = zeros(size(T.time));
T.Perim = zeros(size(T.time));
T.sliver_prev = 0;
T.sliver_i = 0;
T.sliv_fig = 1;

for i=1:length(T.time)
    %clear previous plot
    clf;
    %plot cross-section at current time
    axis equal;
    axis ([-C.g_rad, C.g_rad, -C.g_rad, C.g_rad]);
    set(gca,'units','normalized');
    P.time_str = strcat('Time[s]: ',num2str(T.time(i)));
    xlabel(P.time_str,'FontSize',18);
    T.hole_data(:,3) = T.r(:,i);
    drawcircle(0,0,C.g_rad,'k',1);                      % draw outer shell
    drawholes(T.hole_data,T.num_holes,'w',1);           % draw inner holes
    
    %calculate current port area using binary image processing
    [T.Aport(i), T.Perim(i), T.Sliver, T.sliver_occured] = getHoleArea(C.g_rad);  
    %inches^2, inches, boolean
    
    %pause image if slivers occur
    if (T.sliver_occured == 1 && T.sliver_prev == 0)
        T.sliv_fig = T.sliv_fig+1;
        T.sliver_i = i;
        picture = copyobj(gcf,0);
        figure(bw_plot);
    end
    
    T.sliver_prev = T.sliver_occured;
    
    T.Aport_inches(i) = T.Aport(i);
    T.Aport_inches_act(i) = sum(pi.*T.r(:,i).^2);   % mathematical area wo/ overlap for comparison
    T.Perim_act(i) = sum(2*pi.*T.r(:,i));           % mathematical perimeter wo/ overlap for comparison
    
    %calculate current regression rate (dr)
    T.Aport(i) = T.Aport(i) * (0.0254)^2;  % [m^2]
    T.G = (T.mdotO2(i)/T.Aport(i)); % kg/(m^2s)
    T.dr = F.a*(T.G^F.n)/2.54;  %inches/s
    T.r(:,i+1) = T.r(:,i)+(T.dr*C.time_step);
    disp(T.r)
    
    T.m_burned(i) = C.FG_length*C.p_HDPE*(T.Aport(i)-T.Aport(1)); 
    hold off;
end

hold on;
T.hole_data(:,3) = T.r(:,1);
drawholes(T.hole_data,T.num_holes,'r',0); %draw initial hole size in red

T.fuel_burned = 1+T.sliv_fig;
T.mix_rat = T.fuel_burned + 1;
T.area_port = T.mix_rat + 1;
T.perim_port = T.area_port + 1;

P.p = polyfit(T.time,T.m_burned,2);
P.m_burned_fit = polyval(P.p,T.time);
P.m_burned_fit(P.m_burned_fit<0)=0;

%Plot mass burned over time
% figure(T.fuel_burned);
% plot(T.time, T.m_burned*1000,'r',T.time, P.m_burned_fit*1000,...
%     'r-.','LineWidth',1.2);
% xlabel('time [s]');
% ylabel('fuel burned [g]');
% legend('calculated', 'quadratic fit', 'Location', 'best');

%Plot ideal mixture ratio over time
% figure;
% axis([0,C.burn_time,0,4]);
P.mdot_fuel = (T.m_burned)./T.time;
P.mdot_fuel_fit = (P.m_burned_fit)./T.time;
P.mixture_ratio = T.mdotO2./ P.mdot_fuel;
P.mixture_ratio_fit = T.mdotO2./ P.mdot_fuel_fit;
P.mixture_ratio(isnan(P.mixture_ratio)) = 0;
P.mixture_ratio_fit(isnan(P.mixture_ratio)) = 0;
C.mr_Ave(j) = mean(P.mixture_ratio);
C.surfaceArea(j) = T.Perim_act(1) * C.in2m * C.length;

% plot(T.time,P.mixture_ratio,'b',T.time, P.mixture_ratio_fit,...
%     'b-.','LineWidth',1.2);
% xlabel('time [s]');
% ylabel('mixture ratio [mdotO2/mdotfuel]');
% legend('calculated', 'quadratic fit','Location', 'best');
% P.mixture_ratio;
% T.m_burned;

%Plot total port area over time
% figure(T.area_port);
% T.Aport_inches;
% T.Aport_inches_act;
% 
% P.p2 = polyfit(T.time,T.Aport_inches,2);
% P.Aport_inches_fit = polyval(P.p2,T.time);
% 
% plot(T.time, T.Aport_inches,'k',T.time, P.Aport_inches_fit,'k-.',...
%         T.time, T.Aport_inches_act,'b-.','LineWidth',1.2);
% xlabel('time [s]');
% ylabel('Hole Area [inches^2]');
% legend('With overlap (binary)', 'Quadratic fit(binary)', 'Without overlap',...
%     'Location', 'best');
% 
% P.p3 = polyfit(T.time,T.Perim,2);
% P.Perim_fit = polyval(P.p3,T.time);

%Plot total port perimeter over time
% figure(T.perim_port);
% plot(T.time, T.Perim,'k',T.time,P.Perim_fit,'k-.',T.time, T.Perim_act,...
%         'r-.','LineWidth',1.2);
% xlabel('time [s]');
% ylabel('Hole Perimeter [inches]');
% legend('With overlap (binary)', 'Quadratic fit(binary)', 'Without overlap',...
%     'Location', 'best');

plotfixer;
end

%% Plots of mixture ratio analyses for XOXO_V2 design

E.SA = C.surfaceArea;
E.SA4 = 0.0386;
E.SA5 = 0.04580;
E.SA7 =  0.04580;

E.mr = C.mr_Ave;
E.mr4 = 3.0913;
E.mr5 = 4.1312;
E.mr7 = 3.4609;

%Plot mixture ratio versus surface area
figure;
hold on;
plot(C.surfaceArea, C.mr_Ave, '--');
plot(E.SA4, E.mr4, 'd');
plot(E.SA5, E.mr5, 'o');
plot(E.SA7, E.mr7, 's');
ylabel('MR');
xlabel('Surface Area [m^2]');
legend('Theoretical', 'Test 4',  'Test 5', 'Test 7', 'Location', 'best');
save('mrVsSurfaceAreaData', '-struct', 'E', 'SA', 'SA4', 'SA5', 'SA7', 'mr', 'mr4', 'mr5', 'mr7');
plotfixer
savefig('mrVsSurfaceArea.fig');