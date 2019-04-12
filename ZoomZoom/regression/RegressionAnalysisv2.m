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
burn_time = 3.4;    	% PLACEHOLDER NEED TO BE FIXED
time_step = 0.1;        % PLACEHOLDER NEED TO BE FIXED
pause_time = 0.025;
inches_to_m = 0.0254;
bw_plot = 1

%Get data points from user
%FG_data = xlsread('FuelGrain_Dim7');
%outer_rad = input('Enter Fuel Grain Radius (inches):');
filename = input('Enter .mat file of design: '); 
FG_data = load(filename,'-mat');
outer_rad = 1.005;

hole_data(:,1) = FG_data.x;%FG_data(i,1);
hole_data(:,2) = FG_data.y;%FG_data(i,2);
hole_data(:,3) = FG_data.r;%FG_data(i,3);
num_holes = length(hole_data);

%Start combustion
load('hotplay2.mat');

time = 0:time_step:burn_time;
% time = resize(time,6);
% time_step = time(2)-time(1);

chamP_abs = (101325+(chamP.*6894.7)).*1E-6; %MPa - lab

%get average chamP over firing

thrust(thrust < 5) = 0;
for i=1:length(thrust)
    if thrust(i) == 0
        chamP_abs(i) = 0;
    end
end
chamP_abs(chamP_abs == 0) = [];
% Chamber pressure assumed to be mean value from lab data, can be changed
% for higher accuracy.
P = ones(size(time)).*mean(chamP_abs);  %MPa  

% mdotO2 assumed to be mean value from lab data. 
mdotO2 = ones(size(time))*mean(m_dot_O2);  %kg/s
% mdotO2 = m_dot_O2;

figure(bw_plot);
hold on;

%Draw hole at each time interval with pause time of 0.5 sec
%use equation r=a(G)^n [cm/s] where G=(mdot_O2/Aport) [kg/m^2s]
%coefficients interpolated over lab data for higher accuracy

% possible values of a: 6.3E-3; 1.953e-5; 1.07e-4; 8.3095e-05; 1.7595e-05
% a=4.2908e-05;%6.3E-3;
% n=0.5675; %0.40;
a =  3.5958e-05; % m/s  3.5958e-05 from Sarah
n = 0.6845;


%convert inches to cm
p_HDPE = 950; %kg/m^4
FG_length = 5.350 * (2.54/100); %m

%init arrays
%r(:,1) = hole_data(:,3);
r(:,1) = FG_data.r;  
m_burned = zeros(size(time));
Aport = zeros(size(time));
Aport_inches = zeros(size(time));
Perim = zeros(size(time));
sliver_prev = 0;
sliver_i = 0;
sliv_fig = 1;
%dr = 0;
%rad = outer_rad/2;

for i=1:length(time)
    %clear previous plot
    clf;
    %plot cross-section at current time
    axis equal;
    axis ([-outer_rad, outer_rad, -outer_rad, outer_rad]);
    set(gca,'units','normalized');
    time_str = strcat('Time[s]: ',num2str(time(i)));
    xlabel(time_str,'FontSize',18);
    %set(gca,'Color','k');
    hole_data(:,3) = r(:,i);
    drawcircle(0,0,outer_rad,'k',1);  %draw outer shell
    drawholes(hole_data,num_holes,'w',1); %draw inner holes
    
    %calculate current port area using binary image processing
    [Aport(i), Perim(i), Sliver, sliver_occured] = getHoleArea(outer_rad);  %inches^2, inches, boolean
    
    %pause image if slivers occur
    if (sliver_occured == 1 && sliver_prev == 0)
        sliv_fig = sliv_fig+1;
        sliver_i = i;
        picture = copyobj(gcf,0);
        figure(bw_plot);
    end
    
    sliver_prev = sliver_occured;
    
    Aport_inches(i) = Aport(i);
    Aport_inches_act(i) = sum(pi.*r(:,i).^2);   %mathematical area wo/ overlap for comparison
    Perim_act(i) = sum(2*pi.*r(:,i));           %mathematical perimeter wo/ overlap for comparison
    
    %calculate current regression rate (dr)
    Aport(i) = Aport(i) * (2.54*2.54)/10000;  %m^2
    G = (mdotO2(i)/Aport(i)); % kg/(m^2s)
    dr = a*(G^n)/2.54*100; %(a*(G^n)*(P(i)^b))/2.54;  %inches/s
    r(:,i+1) = r(:,i)+(dr*time_step);
    disp(r)
    
    m_burned(i) = FG_length*p_HDPE*(Aport(i)-Aport(1));  
    pause(pause_time);
    hold off;
end

hold on;
hole_data(:,3) = r(:,1);
drawholes(hole_data,num_holes,'r',0); %draw initial hole size in red

fuel_burned = 1+sliv_fig;
mix_rat = fuel_burned + 1;
area_port = mix_rat + 1;
perim_port = area_port + 1;

p = polyfit(time,m_burned,2);
m_burned_fit = polyval(p,time);
m_burned_fit(m_burned_fit<0)=0;

%Plot mass burned over time
figure(fuel_burned);
plot(time,m_burned*1000,'r',time,m_burned_fit*1000,'r-.','LineWidth',1.2);
xlabel('time [s]');
ylabel('fuel burned [g]');
legend('calculated', 'quadratic fit', 'Location', 'best');

%Plot ideal mixture ratio over time
figure(mix_rat);
axis([0,burn_time,0,4]);
mdot_fuel = (m_burned)./time;
mdot_fuel_fit = (m_burned_fit)./time;
mixture_ratio = mdotO2./ mdot_fuel;
mixture_ratio_fit = mdotO2./ mdot_fuel_fit;
mixture_ratio(isnan(mixture_ratio)) = 0;
mixture_ratio_fit(isnan(mixture_ratio)) = 0;

plot(time,mixture_ratio,'b',time,mixture_ratio_fit,'b-.','LineWidth',1.2);
xlabel('time [s]');
ylabel('mixture ratio [mdotO2/mdotfuel]');
legend('calculated', 'quadratic fit','Location', 'best');
mixture_ratio;
m_burned;

%Plot total port area over time
figure(area_port);
Aport_inches;
Aport_inches_act;

p2 = polyfit(time,Aport_inches,2);
Aport_inches_fit = polyval(p2,time);

plot(time,Aport_inches,'k',time,Aport_inches_fit,'k-.',time,Aport_inches_act,'b-.','LineWidth',1.2);
xlabel('time [s]');
ylabel('Hole Area [inches^2]');
legend('With overlap (binary)', 'Quadratic fit(binary)', 'Without overlap',...
    'Location', 'best');

p3 = polyfit(time,Perim,2);
Perim_fit = polyval(p3,time);

%Plot total port perimeter over time
figure(perim_port);
plot(time,Perim,'k',time,Perim_fit,'k-.',time,Perim_act,'r-.','LineWidth',1.2);
xlabel('time [s]');
ylabel('Hole Perimeter [inches]');
legend('With overlap (binary)', 'Quadratic fit(binary)', 'Without overlap',...
    'Location', 'best');

plotfixer;

% figure(sliver);
% plot(time,Sliver,'r',time,'LineWidth',1.2);
% xlabel('time [s]');
% ylabel('Sliver Area [inches^2]');