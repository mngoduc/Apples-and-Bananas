%% Fuel Grain 
% This script will help design the fuel grain
% We find regression rate and change in grain thickness
clear all; close all;clc;
ourData1 =  'testfire1.mat';
ourData2 =  'testfire2.mat';
ourData3 =  'testfire3.mat';
ourData4 =  'testfire4.mat';
ourData5 =  'testfire5.mat';
ourData6 =  'testfire6.mat';
ourData7 =  'testfire7.mat';
D1 = load(ourData1);
D2 = load(ourData2);
D3 = load(ourData3);
D4 = load(ourData4);
D5 = load(ourData5);
D6 = load(ourData6);
D7 = load(ourData7);

in2m = 0.0254;                      % conversion of inches to meters
C.rho_HDPE = 970;                   % [kg/m^3]  density of HDPE
C.length = 5.375 * in2m;                % [m] length of HDPE

%%
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

D1.A_port_xs = (0.75/2)^2 * pi;
D1.A_port_xs = D1.A_port_xs * in2m^2;

D1.P_port = (0.75) * pi;
D1.P_port = D1.P_port * in2m;
D1.A_surface = D1.P_port * C.length;

D1.r_dot = D1.mdot/(C.rho_HDPE * D1.A_surface); %[kg/s]
D1.d_thick = D1.r_dot * D1.runTime;     % change in wall thickness

D1.mdot_O2 = D1.m_O2/D1.runTime;        % [kg/s]
D1.Go = D1.mdot_O2/(D1.A_port_xs);       % [kg/s/m^2]

%%
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

D2.A_port_xs = (0.75/2)^2 * pi;
D2.A_port_xs = D2.A_port_xs * in2m^2;

D2.P_port = (0.75) * pi;
D2.P_port = D2.P_port * in2m;
D2.A_surface = D2.P_port * C.length;

D2.r_dot = D2.mdot/(C.rho_HDPE *D2.A_surface);
D2.d_thick = D2.r_dot * D2.runTime;         % change in wall thickness

D2.mdot_O2 = D2.m_O2/D2.runTime;
D2.Go = D2.mdot_O2/(D2.A_port_xs);

%%
D3.t1 = D3.time(D3.start_index);
D3.t2 = D3.time(D3.final_index);
D3.dt = D3.time(2) - D3.time(1);
D3.t1i = D3.t1/D3.dt;
D3.t2i = D3.t2/D3.dt;
D3.I = 0;
for i = D3.t1i:D3.t2i
    D3.I = D3.I + D3.thrust(i)*D3.dt;   % [Ns]
end
D3.runTime = D3.t2 - D3.t1;             % [s] time of fire
D3.mdot = D3.mfuel/D3.runTime/1e3;      % [kg/s]

D3.A_port_xs = (0.25/2)^2*pi + 7*(0.182/2)^2*pi + 15*(0.1718/2)^2*pi;
D3.A_port_xs = D3.A_port_xs * in2m^2;

D3.P_port = (0.25)*pi + 7*(0.182)*pi + 15*(0.1718)*pi;
D3.P_port = D3.P_port * in2m;
D3.A_surface = D3.P_port * C.length;

D3.r_dot = D3.mdot/(C.rho_HDPE *D3.A_surface); %[kg/s]
D3.d_thick = D3.r_dot * D3.runTime;     % change in wall thickness

D3.mdot_O2 = D3.m_O2/D3.runTime;        % [kg/s]
D3.Go = D3.mdot_O2/(D3.A_port_xs);       % [kg/s/m^2]

%%
D4.t1 = D4.time(D4.start_index);
D4.t2 = D4.time(D4.final_index);
D4.dt = D4.time(2) - D4.time(1);
D4.t1i = D4.t1/D4.dt;
D4.t2i = D4.t2/D4.dt;
D4.I = 0;
for i = D4.t1i:D4.t2i
    D4.I = D4.I + D4.thrust(i)*D4.dt;   % [Ns]
end
D4.runTime = D4.t2 - D4.t1;             % [s] time of fire
D4.mdot = D4.mfuel/D4.runTime/1e3;      % [kg/s]

D4.A_port_xs = (0.25/2)^2*pi + 4*(0.189/2)^2*pi + 8*(0.1562/2)^2*pi + 4*(0.12/2)^2*pi + 8*(0.1015/2)^2*pi;
D4.A_port_xs = D4.A_port_xs * in2m^2;

D4.P_port = (0.25)*pi + 4*(0.189)*pi + 8*(0.1562)*pi + 4*(0.12)*pi + 8*(0.1015)*pi;
D4.P_port = D4.P_port * in2m;
D4.A_surface = D4.P_port * C.length;

D4.r_dot = D4.mdot/(C.rho_HDPE *D4.A_surface); %[kg/s]
D4.d_thick = D4.r_dot * D4.runTime;     % change in wall thickness

D4.mdot_O2 = D4.m_O2/D4.runTime;        % [kg/s]
D4.Go = D4.mdot_O2/(D4.A_port_xs);       % [kg/s/m^2]

%%
D5.t1 = D5.time(D5.start_index);
D5.t2 = D5.time(D5.final_index);
D5.dt = D5.time(2) - D5.time(1);
D5.t1i = int16(D5.t1./D5.dt);
D5.t2i = int16(D5.t2./D5.dt);
D5.I = 0;
for i = D5.t1i:D5.t2i
    D5.I = D5.I + D5.thrust(i).*D5.dt;   % [Ns]
end
D5.runTime = D5.t2 - D5.t1;             % [s] time of fire
D5.mdot = D5.mfuel/D5.runTime/1e3;      % [kg/s]

D5.A_port_xs = (0.375/2)^2*pi + 4*(0.221/2)^2*pi + 8*(0.1695/2)^2*pi + 4*(0.125/2)^2*pi + 8*(0.136/2)^2*pi;
D5.A_port_xs = D5.A_port_xs * in2m^2;

D5.P_port = (0.375)*pi + 4*(0.221)*pi + 8*(0.1695)*pi + 4*(0.125)*pi + 8*(0.136)*pi;
D5.P_port = D5.P_port * in2m;
D5.A_surface = D5.P_port * C.length;

D5.r_dot = D5.mdot/(C.rho_HDPE *D5.A_surface); %[kg/s]
D5.d_thick = D5.r_dot * D5.runTime;     % change in wall thickness

D5.mdot_O2 = D5.m_O2/D5.runTime;        % [kg/s]
D5.Go = D5.mdot_O2/(D5.A_port_xs);       % [kg/s/m^2]

%%
D6.t1 = D6.time(D6.start_index);
D6.t2 = D6.time(D6.final_index);
D6.dt = D6.time(2) - D6.time(1);
D6.t1i = D6.t1/D6.dt;
D6.t2i = D6.t2/D6.dt;
D6.I = 0;
for i = D6.t1i:D6.t2i
    D6.I = D6.I + D6.thrust(i)*D6.dt;   % [Ns]
end
D6.runTime = D6.t2 - D6.t1;             % [s] time of fire
D6.mdot = D6.mfuel/D6.runTime/1e3;      % [kg/s]

D6.A_port_xs = (0.75/2)^2 * pi;
D6.A_port_xs = D6.A_port_xs * in2m^2;

D6.P_port = (0.75) * pi;
D6.P_port = D6.P_port * in2m;
D6.A_surface = D6.P_port * C.length;

D6.r_dot = D6.mdot/(C.rho_HDPE *D6.A_surface); %[kg/s]
D6.d_thick = D6.r_dot * D6.runTime;     % change in wall thickness

D6.mdot_O2 = D6.m_O2/D6.runTime;        % [kg/s]
D6.Go = D6.mdot_O2/(D6.A_port_xs);       % [kg/s/m^2]

%%
D7.t1 = D7.time(D7.start_index);
D7.t2 = D7.time(D7.final_index);
D7.dt = D7.time(2) - D7.time(1);
D7.t1i = D7.t1/D7.dt;
D7.t2i = D7.t2/D7.dt;
D7.I = 0;
for i = D7.t1i:D7.t2i
    D7.I = D7.I + D7.thrust(i)*D7.dt;   % [Ns]
end
D7.runTime = D7.t2 - D7.t1;             % [s] time of fire
D7.mdot = D7.mfuel/D7.runTime/1e3;      % [kg/s]

D7.A_port_xs = (0.375/2)^2*pi + 4*(0.221/2)^2*pi + 8*(0.1695/2)^2*pi ...
                + 4*(0.125/2)^2*pi + 8*(0.136/2)^2*pi;
D7.A_port_xs = D7.A_port_xs * in2m^2;

D7.P_port = (0.375)*pi + 4*(0.221)*pi + 8*(0.1695)*pi ...
                + 4*(0.125)*pi + 8*(0.136)*pi;
D7.P_port = D7.P_port * in2m;
D7.A_surface = D7.P_port * C.length;

D7.r_dot = D7.mdot/(C.rho_HDPE *D7.A_surface); %[kg/s]
D7.d_thick = D7.r_dot * D7.runTime;     % change in wall thickness

D7.mdot_O2 = D7.m_O2/D7.runTime;        % [kg/s]
D7.Go = D7.mdot_O2/(D7.A_port_xs);       % [kg/s/m^2]

%%

I.n = linspace(0.4,0.8,100);
a = 0;
n = 0;
Rsq = 0;
r_Dot = [D3.r_dot, D4.r_dot, D5.r_dot, D7.r_dot];

for i = I.n
    GoN = [D3.Go^i, D4.Go^i, D5.Go^i, D7.Go^i];
    coef = polyfit(GoN, r_Dot, 1);
    estimData = polyval(coef,GoN);
    R = corrcoef(r_Dot,estimData);
    newRsq = R(1,2).^2;
    if (newRsq - Rsq) > 0.001
        a = coef(1);
        n = i;
        Rsq = newRsq;
    end
end

figure;
hold on;
plot(GoN, r_Dot, '*');
plot(GoN, estimData);
xlabel("G_o^N");
ylabel("R_{Dot}");
legend("Experimental Data", "Linear Fit", "location", "southeast");
plotfixer;

filename = 'regression';
save(filename)