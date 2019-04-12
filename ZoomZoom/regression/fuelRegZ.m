%% Fuel Grain Function
% This script will help design the fuel grain
% We find regression rate and change in grain thickness
function [ a_1, n ] = fuelRegZ( ourData )
clear all; close all;clc;
ourData =  'testfire4.mat';
isabelData = 'RocketTestFire.mat';
D1 = load(ourData);
D2 = load(isabelData);

C.rho_HDPE = 970;                   % [kg/m^3]  density of HDPE
C.r = 0.75*2.54/100;%0.009525;      % [m] radius of hole
C.L = 5.375*2.54/100;               % [m] length of fuel grain (5.375 in)
C.A_surface = 2*pi*C.r*C.L;         % [m^2]<---- place holder surface area of grain
C.A_port_xs = pi*C.r^2;               % cross sectional area of ports


D1.t1 = D1.time(D1.start_index);
D1.t2 = D1.time(D1.final_index);
D1.dt = D1.time(2) - D1.time(1);
D1.t1i = D1.t1/D1.dt;
D1.t2i = D1.t2/D1.dt;
D1.I = 0;
for i = D1.t1i:D1.t2i
    D1.I = D1.I + D1.thrust(i)*D1.dt;   % [Ns]
end
D1.runTime = D1.t2 - D1.t1;             % [s] time of fire
D1.mdot = D1.mfuel/D1.runTime/1e3;      % [kg/s]

D1.r_dot = D1.mdot/(C.rho_HDPE *C.A_surface); %[kg/s]
D1.d_thick = D1.r_dot * D1.runTime;     % change in wall thickness

D1.A_port_xs = 0.125^2*pi + 7*(0.185/2)^2*pi + 15*(0.1718/2)^2*pi;
D1.A_port_xs = D1.A_port_xs *0.0254^2;
D1.mdot_O2 = D1.m_O2/D1.runTime;        % [kg/s]
D1.Go = D1.mdot_O2/(D1.A_port_xs);       % [kg/s/m^2]

D2.t1 = D2.time(D2.start_index);
D2.t2 = D2.time(D2.final_index);
D2.dt = D2.time(2) - D2.time(1);
D2.t1i = D2.t1/D2.dt;
D2.t2i = D2.t2/D2.dt;

D2.I = 0;
for i = D2.t1i:D2.t2i
    D2.I = D2.I + D2.thrust(i)*D2.dt; % [Ns]
end

D2.runTime = D2.t2 - D2.t1;             % [seconds] time of fire
D2.mdot = D2.mfuel/D2.runTime/1e3;      % [kg/s]

D2.r_dot = D2.mdot/(C.rho_HDPE *C.A_surface);
D2.d_thick = D2.r_dot * D2.runTime;         % change in wall thickness

D2.mdot_O2 = D2.m_O2/D2.runTime;
D2.Go = D2.mdot_O2/(C.A_port_xs);


I.n = 0.38;
I.dn= 1e-4;

I.a_1 = D1.r_dot/D1.Go^I.n;
I.a_2 = D2.r_dot/D2.Go^I.n;

while abs(I.a_1 - I.a_2) > 2.5e-5
    I.n = I.n + I.dn;
    
    I.a_1 = D1.r_dot/D1.Go^I.n;
    I.a_2 = D2.r_dot/D2.Go^I.n;
    
    if I.n > 2
        displ(['n got too big, n = ', n])
        break
    end
end

a_1 = I.a_1;
n = I.n;
end