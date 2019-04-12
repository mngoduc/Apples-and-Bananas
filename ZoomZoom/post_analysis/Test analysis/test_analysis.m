
%% Test Analysis
% ME 140 - Reginald Mitchell - Spring 2018
% Projects 7 - The Hybrid Rocket Motor
% Alex Casa, Ian Gunady, Minh Ngo Duc, Thomas Trzpit, Zane Zook

% analyze the test fires
% testfire1 & testfire2 are from Isabel and Sarah
% ===============================================
% Important variables from Lab 3 data sets (testfire1 and testfire2)
% A_throat  = 0.0001824 [m^2]
% P_choked_min = choked orifice pressure (kPa gage)
% chamP = chamber pressure [kPa gage]
% resP = O2 resevoir pressure (kPa gage)
% thrust = thrust [N]
% temp = temp [K]
% mfuel = mass of consumed fuel [g]
% m_O2 = mass O2 consumed [kg]
% m_total = mass O2 + mass fuel consumed [kg]
% ==================================================
% 5/26/18 3:30pm testfire3 our own fuel grain, 5/32'' and 1/1'' hole diameter. 5.375'' to
% tal length, with upstream mixing chamber
% 5/31/18 2:00pm testfire4 our own fuel grain, X pattern, cone mixing
% chamber, 5.375'' total length
% 6/3/18 1:45pm testfire5 our own fuel grain, X pattern, cone mixing
% chamber, 5.375'' total length, delayed ignition
% 6/4/18 11:30am testfire6, TA fuel grain, our nozzzle
% 6/5/18 1:15pm testfire7, our own fuel grain, X pattern, cone mixing
% chamber, 5.375'' total length, our own converging (5th order polynomial)
% and diverging (conical) nozzle
clear all; close all; clc;
A = load('testfire1.mat');
B = load('testfire2.mat');
C = load('testfire3.mat');
D = load('testfire4.mat');
E = load('testfire5.mat');
F = load('testfire6.mat');
G = load('testfire7.mat');

% foul? 0 if did not foul, 1 if did foul
A.foul = 0;
B.foul = 0;
C.foul = 0;
D.foul = 1;
E.foul = 0;
F.foul = 0;
G.foul = 0;
set(0,'DefaultAxesColor','none')
% plot thrust vs time
plot(A.time,A.thrust)
hold on
plot(B.time,B.thrust)
plot(C.time,C.thrust)
plot(D.time,D.thrust)
plot(E.time,E.thrust)
plot(F.time,F.thrust)
plot(G.time,G.thrust)
xlabel('Time [s]'); ylabel('Thrust [N]')
legend('TA data 1','TA data 2','test 1','test 2','test 3','test 4','test 5',...
    'Location', 'best');
plotfixer;

% plot pressures vs time

% calculate total and specific impulse
g = 9.81;   % [m/s^2]
dt = 0.001; % [s]
Best.impulse = 891;

A.mfuel = A.mfuel/1000; % [kg]
B.mfuel = B.mfuel/1000;
C.mfuel = C.mfuel/1000;
D.mfuel = D.mfuel/1000;
E.mfuel = E.mfuel/1000;
F.mfuel = F.mfuel/1000;
G.mfuel = G.mfuel/1000;

Best.impulse_sp = 218;

% calculate mixture ratio
A.mr = A.m_O2/A.mfuel;
B.mr = B.m_O2/B.mfuel;
C.mr = C.m_O2/C.mfuel;
D.mr = D.m_O2/D.mfuel;
E.mr = E.m_O2/E.mfuel;
F.mr = F.m_O2/F.mfuel;
G.mr = G.m_O2/G.mfuel;

ThrustAll.A = A.thrust;
ThrustAll.B = B.thrust;
ThrustAll.C = C.thrust;
ThrustAll.D = C.thrust;
ThrustAll.E = E.thrust;
ThrustAll.F = F.thrust;
ThrustAll.G = G.thrust;
ThrustAll.time = A.time;
save('ThrustAll', '-struct', 'ThrustAll', 'time', 'A', 'B', 'C', 'D', 'E', 'F', 'G')

% Charactersitic area ratio - cross sectional to surface area ratio
% test1 and test2 (assuming a single hole with 3/4'' diameter)
A.r = (19.05/2)*1e-3;       % [m] radius of circle
A.L = 136.525*1e-3;         % [m] length of grain (5.375 in)
A.ax = pi*A.r^2;            % [m^2] initial total cross sectional area of port (~O2)
A.as = 2*pi*A.r;            % [m^2] initial total surface area of port (~fuel)
A.A_rat = A.ax/A.as;        % characteristic area ratio (cross sectional to surface area ratio)

% test3
C.rs = (3.96875/2)*1e-3;    % [m] radius of smaller circles
C.rl = (6.35/2)*1e-3;       % [m] radius of the larger circle in the middle
C.L = 136.525*1e-3;         % [m] length of grain (5.375 in)
C.ax = pi*C.rl^2 + 22*pi*C.rs^2;    % [m^2] initial total cross sectional area of ports (~O2)
C.as = 2*pi*C.rl*C.L + 22*2*pi*C.rs;% [m^2] initial total surface area of ports (~fuel)
C.A_rat = C.ax/C.as;        % characteristic area ratio (cross sectional to surface area ratio)
%*85
% Table of results
TestName = ['TA data 1';'TA data 2';'test 1   ';'test 2   ';'test 3   ';'test 4   ';'test 5   '];
mr = [A.mr;B.mr;C.mr;D.mr;E.mr;F.mr;G.mr];
Impulse = [A.impulse;B.impulse;C.impulse;D.impulse;E.impulse;F.impulse;G.impulse];
SpecificImpulse = [A.Isp;B.Isp;C.Isp;D.Isp;E.Isp;F.Isp;G.Isp];
Foul = [A.foul;B.foul;C.foul;D.foul;E.foul;F.foul;G.foul];

T = table(TestName,mr,Impulse,SpecificImpulse,Foul);
disp(T)

%% Final test analysis
Patm = 101325;      % [Pa]

N.dt = 0.0166624;   % [m] nozzle throat diameter (0.66 inch)
N.rt = N.dt/2;      % [m] nozzle throat radius
N.At = pi*N.rt^2;   % [m^2] nozzle throat area 
N.epsilon = 1.81;   % expansion ratio <- DOUBLE CHECK THIS !!!!!!
N.Ae = N.epsilon*N.At;  % [m^2] exit area

G.dt = G.time(G.final_index) - G.time(G.start_index);   % [s] time of test
G.thrust_avg = G.impulse/G.dt;      % [N]
G.mdot = G.m_total/G.dt;            % [kg/s]
G.v_exit_avg = G.thrust_avg/G.mdot; % [m/s] thrust = ve  + pressure (assume pressure terms cancel)
G.Po = G.chamP(G.start_index:G.final_index);    % [kPa] gage
G.Po = G.Po*1000 + Patm;            % [Pa]
G.Po_avg = mean(G.Po);              % [Pa]

G.Me = 2.83;    % Mach at exit, based on nozzle expansion ratio, assumes isentropic
% M = Ve/c
G.ce = G.v_exit_avg/G.Me;  % [m/s] speed of sound at exit
% c = krt
% G.Te = G.ce/G.k*R;

% %% Pressure test B
% figure
% plot(B.time,B.chamP)
% hold on
% plot(B.time,B.P_choked_min)
% legend('chamber pressure','resevoir pressure'); title('test B pressures')
% xlabel('time'); ylabel('pressure')
% plotfixer
% 
% %% Pressure test C
% figure
% plot(C.time,C.chamP)
% hold on
% plot(C.time,C.P_choked_min)
% legend('chamber pressure','resevoir pressure'); title('test C pressures')
% xlabel('time'); ylabel('pressure')
% plotfixer
% 
%% Pressure test D
figure
plot(D.time,D.chamP)
hold on
plot(D.time,D.P_choked_min)
legend('chamber pressure','resevoir pressure'); title('test D pressures')
xlabel('Time [s]'); ylabel('Pressure [kPa] gage')
plotfixer

%% total impulse vs mr
Imr.mr = [A.mr B.mr C.mr D.mr E.mr F.mr G.mr];
Imr.impulse = [A.impulse B.impulse C.impulse D.impulse E.impulse F.impulse G.impulse];
save('Impulse_vs_mr', '-struct', 'Imr', 'mr', 'impulse')
plot(Imr.mr,Imr.impulse,'o')
xlabel('Mixture Ratio')
ylabel('Total Impulse [Ns]')
plotfixer

%% Pressure test E
E.dt = E.time(E.final_index) - E.time(E.start_index);   % [s] time of test
E.thrust_avg = E.impulse/E.dt;      % [N]
E.mdot = E.m_total/E.dt;            % [kg/s]
E.v_exit_avg = E.thrust_avg/E.mdot; % [m/s] thrust = ve  + pressure (assume pressure terms cancel)
E.Po = E.chamP(E.start_index:E.final_index);    % [kPa] gage
E.Po = E.Po*1000 + Patm;            % [Pa]
E.Po_avg = mean(E.Po);              % [Pa]
figure
plot(E.time,E.chamP)
hold on
plot(E.time,E.P_choked_min)
legend('chamber pressure','resevoir pressure'); title('test E pressures')
xlabel('Time [s]'); ylabel('Pressure [kPa] gage')
plotfixer
save('Pressure5', '-struct', 'E', 'time', 'P_choked_min', 'chamP')
%% Pressure test F
figure
plot(F.time,F.chamP)
hold on
plot(F.time,F.P_choked_min)
legend('chamber pressure','resevoir pressure'); title('test F pressures')
xlabel('Time [s]'); ylabel('Pressure [kPa] gage')
plotfixer
%% Pressure test G
figure
plot(G.time,G.chamP)
hold on
plot(G.time,G.P_choked_min)
legend('chamber pressure','resevoir pressure'); title('test G pressures')
xlabel('Time [s]'); ylabel('Pressure [kPa] gage')
plotfixer
save('Pressure7', '-struct', 'G', 'time', 'P_choked_min', 'chamP')
% %% improvement over test fires
% figure
% plot(Impulse)
% xlabel('Test Fire #')
% ylabel('Total Impulse [Ns]')
% plotfixer