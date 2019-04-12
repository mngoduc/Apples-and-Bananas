function m_fuel = get_fuel_burned(a,n,b,filename,labdata)
% use equation r=a(G)^n*P^b [cm/s] where G=(mdot_O2/Aport) [kg/m^2s]
% filename is xlsx file with fuel grain geometry
% labdata is mat file for lab data that is being used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%init constants
burn_time = 3.4;
time_step = 0.1;
bw_plot = 1;

%Get data points from user
FG_data = xlsread(filename);
outer_rad = 1.995/2;

s = size(FG_data);
hole_data = zeros(size(FG_data));
num_holes = s(1);

%Create hole data
for i=1:num_holes
    hole_data(i,1) = FG_data(i,1);
    hole_data(i,2) = FG_data(i,2);
    hole_data(i,3) = FG_data(i,3);
end

%Start combustion
load(labdata);

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
P = ones(size(time)).*mean(chamP_abs);  %MPa  - assumed

%P_data = ones(size(time))*mean(chamP)*0.006894759086775369; %MPa - act (constant)
% P = abs(chamP*0.006894759086775369);
% P=resize(P,6);
mdotO2 = ones(size(time))*mean(m_dot_O2);  %kg/s
% mdotO2 = m_dot_O2;
% mdotO2 = resize(mdotO2,6);

figure(bw_plot);
%Draw hole at each time interval with pause time of 0.5 sec

%use equation r=a(G)^n*P^b [cm/s] where G=(mdot_O2/Aport) [kg/m^2s]
%regression rate constants 

%convert inches to cm
p_HDPE = 950; %kg/m^3
FG_length = 5.350 * (2.54/100); %m

%init arrays
r(:,1) = hole_data(:,3);
m_burned = zeros(size(time));
Aport = zeros(size(time));
Aport_inches = zeros(size(time));
Perim = zeros(size(time));

for i=1:length(time)
    %clear previous plot
    clf;
    %plot cross-section at current time
    axis equal;
    axis ([-outer_rad outer_rad -outer_rad outer_rad]);
    set(gca,'units','normalized');
    time_str = strcat('Time[s]: ',num2str(time(i)));
    xlabel(time_str,'FontSize',18);
    hole_data(:,3) = r(:,i);
    drawcircle(0,0,outer_rad,'k',1);  %draw outer shell
    drawholes(hole_data,num_holes,'w',1); %draw inner holes

    %calculate current port area using BW image processing
    [Aport(i), Perim(i), Sliver, sliver_occured] = getHoleArea(outer_rad);  %inches^2, inches
    
    %calculate current regression rate (dr)
    Aport(i) = Aport(i) * (2.54*2.54)/10000;  %m^2
    G = (mdotO2(i)/Aport(i)); % kg/(m^2s)
    dr = (a*(G^n)*(P(i)^b))/2.54;  %inches/s
    r(:,i+1) = r(:,i)+(dr*time_step);
    
    m_burned(i) = FG_length*p_HDPE*(Aport(i)-Aport(1));  
    hold off;
end

p = polyfit(time,m_burned,2);
m_burned_fit = polyval(p,time);
m_burned_fit(m_burned_fit<0)=0;

m_fuel = m_burned_fit(end)*1000;
end

