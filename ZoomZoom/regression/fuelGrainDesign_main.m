%% Fuel Grain
% This script will help design the fuel grain
% We find regression rate and change in grain thickness
clear all; close all;
sarahData =  'testfire1.mat';
isabelData = 'RocketTestFire.mat';
load(sarahData);
plot(time(start_index:final_index),thrust(start_index:final_index));
xlabel('time'); ylabel('thrust'); plotfixer();

rho_HDPE = 970;             % [kg/m^3]  density of HDPE

t1 = time(start_index); 
t2 = time(final_index);

dt = time(2) - time(1); 
t1i = t1/dt; t2i = t2/dt;
impulse = 0;
for i = t1i:t2i
    impulse = impulse + thrust(i)*dt; % [Ns]
end
deltaT = t2 - t1;           % [seconds] time of fire
mdot = mfuel/deltaT/1000;   % [kg/s]
r = 0.75*2.54/100;%0.009525;               % [m] radius of hole
L = 0.136271;               % [m] length of fuel grain (5.375 in)
A_surface = 2*pi*r*L;  % <---- place holder surface area of grain
% mdot = r*rho*A
rdot = mdot/(rho_HDPE*A_surface);   % regression rate
delta_thickness = rdot*deltaT; % change in wall thickness

mdot_O2 = m_O2/deltaT;
A_ports = pi*r^2;
n = 0.4; %0.55;
Go = mdot_O2/A_ports;

Go_n = Go^n;
a_exp = rdot/Go_n;

disp('finished running')

%% large loop of 