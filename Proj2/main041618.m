close all; clear all ; clc; 

%% Setting Up 
in2m = 0.0254;
C2K = 273.15; 
lbf2N = 4.44822; 
P_atm = 101325;             % [Pa]
R = 287;                    % [J/kg K]

G1.A1 = 27.3 * in2m^2;              
G1.A2 = 6.4 * in2m^2;
G1.A3 = 9.0 * in2m^2; 
G1.A4 = 7.2 * in2m^2;
G1.A5 = 4.7 * in2m^2; 
G1.A8 = 3.87 * in2m^2;

G2.RF_T2 = 0.68;           % Cross Flow RF
G2.RF_T3 = 0.68;           % Cross Flow RF 
G2.RF_T4 = 0.68;           % Cross Flow RF 
G2.RF_T5 = 0.86;           % Axial Flow RF 
G2.RF_T8 = 0.68;           % Cross Flow RF 

G3.LHV = 42.8e6;        % J/kg

% Importing Data & Reassigning the columns to corresponding names

filename = 'Data.xlsx';
data = xlsread(filename);

D.RPM = data(:,1)';
D.T2m = data(:,2)' + C2K;                               % [K]
D.T3m = data(:,3)' + C2K;                               % [K]
D.T4m = data(:,4)' + C2K;                               % [K]
D.T5m = data(:,5)' + C2K;                               % [K]
D.T8m = data(:,6)' + C2K;                               % [K]
D.DP2 = (data(:,7)')*1e3; % difference b/w static and stagnation, [Pa]
D.PT3 = (data(:,8)')*1e3 + P_atm;                       % [Pa]
D.P4 = (data(:,9)')*1e3 + P_atm;                        % [Pa]
D.PT5 = (data(:,10)')*1e3 + P_atm;                      % [Pa]
D.PT8 = (data(:,11)')*1e3+ P_atm;                       % [Pa]
D.FuelFlow = data(:,12)';                               % [kg/s]
D.Thrust = data(:,13)'*lbf2N;                           % [N]


%% Problem 2 

% State 2
P2.Po2 = P_atm*ones(size(D.DP2));   % [Pa]
P2.P2 = P2.Po2 - D.DP2;             % [Pa]
[~, P2.k2] = specHeatAir(D.T2m); 
P2.M2 = ((2./(P2.k2-1)).*((P2.Po2./P2.P2).^((P2.k2-1)./P2.k2)-1)).^(1/2);

P2.T2 = Tm2T(D.T2m, G2.RF_T2, P2.M2, P2.k2);        % [K]
P2.To2 = T2To(P2.T2, P2.M2, P2.k2);                 % [K]

P2.V2 = P2.M2.*sqrt(P2.k2*R.*P2.T2);
P2.ro2 = P2.Po2./(R.*P2.To2);
P2.r2 = P2.ro2.*(P2.P2./P2.Po2).*(P2.To2./P2.T2);
P2.mdot_air = P2.r2.*P2.V2.*G1.A2;                  % [kg/s]
P2.MFP2 = (P2.mdot_air.*sqrt(R*P2.To2)./(G1.A2.*P2.Po2));

% State 1 
P2.Po1 = P2.Po2;                                    % [Pa]
P2.To1 = P2.To2;                                    % [K]

P2.MFP1 = P2.MFP2*G1.A2/G1.A1;
P2.k1 = P2.k2;
for i = 1:length(P2.k1)
    P2.M1(i) = MFP2M(P2.k1(i),P2.MFP1(i));
end

P2.Tratio1 = isenRatioT(P2.k1, P2.M1);
% using isentropic ratios, returns To/T
P2.T1 = P2.To1./P2.Tratio1;
P2.V1 = P2.M1.*sqrt(P2.k1.*R.*P2.T1);

% State 3
P2.Po3 = D.PT3;
[~, P2.k3] = specHeatAir(D.T3m); 
P2.To3 = D.T3m; % estimated value; first guess

for j = 1:5
    P2.MFP3 = (P2.mdot_air.*sqrt(R*P2.To3)./(G1.A3.*P2.Po3));   % estimated MFP
    for i = 1:length(P2.MFP3)
        P2.M3(i) = MFP2M(P2.k3(i),P2.MFP3(i));   % estimated value
    end
    P2.T3 = D.T3m./(1+G2.RF_T3.*(P2.k3-1)./2.*P2.M3.^2);
    P2.To3 = P2.T3.*(1+(P2.k3-1)./2.*P2.M3.^2);
end

P2.V3 = P2.M3.*sqrt(P2.k3.*R.*P2.T3);

% State 4
P2.P4 = D.P4;
P2.QDotComb = G3.LHV * D.FuelFlow; %[J/s]
P2.mDotTot = D.FuelFlow + P2.mdot_air;

P2.To4 = D.T4m;
[~, P2.k4] = specHeatAir(D.T4m);

for j = 1:5 
    for i = 1:length(P2.To4)
        [P2.MFP4(i),P2.Po4(i),P2.M4(i)] = ToP2MFP(P2.mDotTot(i),R,...
                                    P2.To4(i),G1.A4,P2.k4(i),P2.P4(i));
    end
    P2.T4 = D.T4m./(1+G2.RF_T4.*(P2.k4-1)./2.*P2.M4.^2);
    P2.To4 = P2.T4.*isenRatioT(P2.k4,P2.M4);
    P2.V4 = P2.M4.*sqrt(P2.k4.*R.*P2.T4);
end

% State 5
P2.Po5 = D.PT5;
P2.To5 = D.T5m;  % initilizing to an estimated value

[~, P2.k5] = specHeatAir(D.T5m);
for j = 1:5
    P2.MFP5 = (P2.mDotTot.*sqrt(R*P2.To5)./(G1.A5.*P2.Po5));
    for i = 1:length(P2.k5)
        P2.M5(i) = MFP2M(P2.k5(i),P2.MFP5(i));
    end
    % using isentropic ratios, returns To/T
    P2.T5 = D.T5m./(1+G2.RF_T5.*(P2.k5-1)./2.*P2.M5.^2);
    P2.To5 = P2.T5.*isenRatioT(P2.k5,P2.M5);
end
P2.V5 = P2.M5.*sqrt(P2.k5.*R.*P2.T5);

% State 8
P2.P8 = P2.Po2;         % P8 is the same as Patm because subsonic
P2.Po8 = D.PT8;         
[~, P2.k8] = specHeatAir(D.T8m);
P2.M8 = sqrt(2./(P2.k8-1).*((P2.Po8./P2.P8).^((P2.k8-1)./P2.k8)-1));
P2.T8 = D.T8m./(1+G2.RF_T8.*(P2.k8-1)./2.*P2.M8);
P2.To8 = P2.T8.*isenRatioT(P2.k8,P2.M8);

P2.V8 = (2*dh(P2.To8, P2.T8)).^(1/2);
P2.F_T = P2.mDotTot .* P2.V8; 

%% Plotting Portion of Problem 2

Plot2.Po_arr = [P2.Po1; P2.Po2; P2.Po3; P2.Po4; P2.Po5; P2.Po8];
Plot2.To_arr = [P2.To1; P2.To2; P2.To3; P2.To4; P2.To5; P2.To8];
Plot2.M_arr = [P2.M1; P2.M2; P2.M3; P2.M4; P2.M5; P2.M8];
Plot2.V_arr = [P2.V1; P2.V2; P2.V3; P2.V4; P2.V5; P2.V8];
Plot2.leg = ['Station 1';'Station 2';'Station 3';'Station 4';...
            'Station 5'; 'Station 8'];

% Stagnation Temp Plot
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
for i = 1:6
    Plot2.h(i) = plot(D.RPM/(1e3), Plot2.To_arr(i,:),'o--','LineWidth',2,'MarkerSize',5);
end
legend(Plot2.h,Plot2.leg, 'Location', 'best');
xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
ylabel('Stagnation Temperature, T_o [K]','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);

% Stagnation Pressure Plot 
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
for i = 1:6
    Plot2.h(i) = plot(D.RPM/(1e3), Plot2.Po_arr(i,:)/(1e3),'o--','LineWidth',2,'MarkerSize',5);
end
legend(Plot2.h,Plot2.leg, 'Location', 'best');
xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
ylabel('Stagnation Pressure, P_o [kPa]','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);

% Mach Number Plot 
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
for i = 1:6
    Plot2.h(i) = plot(D.RPM/(1e3), Plot2.M_arr(i,:),'o--','LineWidth',2,'MarkerSize',5);
end
legend(Plot2.h,Plot2.leg, 'Location', 'best');
xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
ylabel('Mach Number, M ','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);

% Velocity Plot 
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
for i = 1:6
    Plot2.h(i) = plot(D.RPM/(1e3), Plot2.V_arr(i,:),'o--','LineWidth',2,'MarkerSize',5);
end
legend(Plot2.h,Plot2.leg, 'Location', 'best');
xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
ylabel('Velocity, V [m/s]','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);

% Air Fuel Ratio:
P2.AF = P2.mdot_air./D.FuelFlow;

% air flow rate, fuel flow rate, 
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(D.RPM/(1e3),P2.mdot_air*1e2, 'o--','LineWidth',2,'MarkerSize',5);
hold on;
plot(D.RPM/(1e3),D.FuelFlow*1e4, 'o--','LineWidth',2,'MarkerSize',5);
plot(D.RPM/(1e3),P2.AF,'o--','LineWidth',2,'MarkerSize',5);
legend('100 x Air Flow Rate', '10000 x Fuel Flow Rate', 'Air-Fuel Ratio', 'Location', 'best')
xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
ylabel('Flow Rates [kg/s] and Air-Fuel Ratio','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);

% calculated thrust and recorded thrust 
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(D.RPM/(1e3),D.Thrust,'o--','LineWidth',2,'MarkerSize',5);
hold on;
plot(D.RPM/(1e3),P2.F_T, 'o--','LineWidth',2,'MarkerSize',5);
legend('Measured Thrust','Calculated Thrust','Location', 'best')
xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
ylabel('Thrust Force, F_T [N]','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);

%% Problem 3

% specific thrust 
P3.ST = P2.F_T ./ P2.mdot_air; 
% [N/kg/s] Specific Thrust, ST, Thrust per unit of air mass flow through engine

% thrust specific fuel consumption
P3.TSFC = D.FuelFlow ./ P2.F_T;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot

yyaxis left
plot(D.RPM/(1e3),P3.ST,'o--','LineWidth',2,'MarkerSize',5);
ylabel('Specific Thrust [N/kg/s]','FontSize',18,'FontWeight','bold');

yyaxis right
plot(D.RPM/(1e3),P3.TSFC*(1e3),'o--','LineWidth',2,'MarkerSize',5);
yyaxis right
ylabel('Thrust Specific Fuel Consumption [mN/kg/s]','FontSize',18,'FontWeight','bold');

xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);
 
% thermal efficiency  = net work / heat input;
% heat input = LHV * m_dot_fuel = P2.QDotComb
% P2.QDotComb = G3.LHV * D.FuelFlow; %[J/s]
% W_net kinetic energy 
P3.W_net = 0.5*P2.mDotTot.*P2.V8.^2;
P3.Eff_therm = P3.W_net./P2.QDotComb; 

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(D.RPM/(1e3),P3.Eff_therm,'o--','LineWidth',2,'MarkerSize',5);
xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
ylabel('Thermal Efficiency, \eta_{therm}','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);

%% Problem 4

P4.W_comp = P2.mdot_air .* dh(P2.To3, P2.To2);
P4.W_turb = -1*P2.mDotTot .* dh(P2.To5, P2.To4);

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(D.RPM/(1e3),P4.W_comp/(1e3),'o--','LineWidth',2,'MarkerSize',5);
hold on; 
plot(D.RPM/(1e3),P4.W_turb/(1e3),'o--','LineWidth',2,'MarkerSize',5);
xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
ylabel('Compressor and Turbine Power, [kW]','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);
legend('Compressor Power', 'Turbine Power', 'Location', 'best');

% "plot of compressor, turbine, and nozzle adiabatic efficiencies, 
% combustor stagnation pressure loss, and 'apparent combustion efficiency' 
% vs spool speed"

% compressor adiabatic efficiencies 
for i = 1:length(P4.W_comp)
    P4.To3s(i) = presRatio2Temp(P2.Po2(i), P2.Po3(i), P2.To2(i));
end
P4.W_comp_s = P2.mdot_air .* dh(P4.To3s, P2.To2);
P4.eta_comp = P4.W_comp_s./P4.W_comp;

% turbine adiabatic efficiencies 
for i = 1:length(P4.W_comp)
    P4.To5s(i) = presRatio2TempInit(P2.Po5(i), P2.Po4(i), P2.To4(i));
end
P4.W_turb_s = -1*P2.mDotTot .* dh(P4.To5s, P2.To4);
P4.eta_turb = P4.W_turb./P4.W_turb_s;

% nozzle adiabatic efficiencies 
P4.T8s = P2.To5.*(P2.P8/P2.Po5).^((P2.k8-1)./P2.k8);
P4.M8s = sqrt((2./(P2.k8 - 1)).*((P2.Po5./P2.P8).^((P2.k8 - 1)./P2.k8)-1));
P4.V8s = P4.M8s.*sqrt(P2.k8.*R.*P4.T8s);
% the way reggie suggested 
P4.V8s_alt = (-2*dh(P2.To5, P2.To8)).^(1/2);

% nozzle efficiency is: 
P4.eta_noz = (P2.V8.^2)./(P4.V8s.^2);

% combustor stagnation pressure loss
P4.CombPressLoss = P2.Po4 ./ P2.Po3;

% apparent combustion efficiency vs spool speed 
P4.Q_comb = P2.mDotTot.* dh(P2.To4, P2.To3);
P4.eta_comb = (P2.mDotTot.* dh(P2.To4, P2.To3))./P2.QDotComb;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(D.RPM/(1e3),P4.CombPressLoss,'o--','LineWidth',2,'MarkerSize',5);
hold on; 
plot(D.RPM/(1e3),P4.eta_comp,'o--','LineWidth',2,'MarkerSize',5);
plot(D.RPM/(1e3),P4.eta_turb,'o--','LineWidth',2,'MarkerSize',5);
plot(D.RPM/(1e3),P4.eta_noz,'o--','LineWidth',2,'MarkerSize',5);
plot(D.RPM/(1e3),P4.eta_comb,'o--','LineWidth',2,'MarkerSize',5);
xlabel('Spool Speed, \omega [kRPM]','FontSize',18,'FontWeight','bold');
ylabel('Efficiency, Pressure Ratio','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);
legend('Apparent Combustion Efficiency', 'Compressor Efficiency',...
    'Turbine Efficiency', 'Nozzle Efficiency', 'Combustor Efficiency',...
    'Location', 'best');
