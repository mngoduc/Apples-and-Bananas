%% Main Script
% ME 140 - Reginald Mitchell - Spring 2018
% Projects 7 - The Hybrid Rocket Motor
% Alex Casa, Ian Gunady, Minh Ngo Duc, Thomas Trzpit, Zane Zook

% x is mol fraction
% y is mass fraction
%% Constants
clear all; close all; clc; 
tic

% C struct - constants
C.P_atm = 100000;           % [Pa] atmospheric pressure, 1 bar
C.P_chamber = 9.7166*100000;    % [Pa] combustion chamber pressure, 68 bar
C.To = 298;                 % [K] standard temperature

% Molar masses
C.M_CO2 = 44.01e-3;         % [kg/mol]
C.M_O2 = 32e-3;             % [kg/mol]
C.M_H2O = 18e-3;            % [kg/mol]
C.M_C2H4 = 28.05e-3;        % [kg/mol]
C.M_C = 12.0107e-3;         % [kg/mol]
C.M_H2 = 2.01588e-3;        % [kg/mol]
C.M_CH2 = C.M_C + C.M_H2;   % [kg/mol]


%% Problem 1 State 1
P1.HHV_HDPE = 46.5e6; %J/kg

P1.hForm_CO2_bar = -393520;                     % [J/mol]
P1.hForm_CO2 = P1.hForm_CO2_bar / C.M_CO2;      % [J/kg]

P1.hForm_O2_bar = 0;                            % [J/mol]
P1.hForm_O2 = 0 / C.M_O2;                       % [J/kg]

P1.hForm_H2O_bar = -285830;                     % [J/mol]
P1.hForm_H2O = P1.hForm_H2O_bar / C.M_H2O;      % [J/kg]

P1.hForm_C2H4_bar = 52280;                      % [J/mol]
P1.hForm_C2H4 = P1.hForm_C2H4_bar / C.M_C2H4;   % [J/kg]

P1.hForm_H2O_norm = C.M_H2O / (C.M_C + C.M_H2) * P1.hForm_H2O;
P1.hForm_CO2_norm = C.M_CO2 / (C.M_C + C.M_H2) * P1.hForm_CO2;
P1.hForm_O2_norm = C.M_O2 / (C.M_C + C.M_H2) * P1.hForm_O2;
P1.hForm_HDPE = (P1.hForm_H2O_norm + P1.hForm_CO2_norm - ...
    P1.hForm_O2_norm) + P1.HHV_HDPE;

P1.dh_shift = P1.hForm_HDPE - P1.hForm_C2H4;    % [J/kg]

%% Problem 1 Combustion
% C2H4 + 3*O2 -> products
P1.gas = IdealGasMix('me140_species.xml');
set(P1.gas,'T',C.To,'P',C.P_chamber);

% Get the Molar Masses from the gas object
MM = molarMasses(P1.gas);

% Create the composition vectors
P1.nsp = nSpecies(P1.gas);

% Finding mass of products
P1.HDPE_mass = 0.085; %kg
P1.HDPE_moles = P1.HDPE_mass / C.M_CH2;
P1.CO2_moles = P1.HDPE_moles;
P1.total_O_moles = P1.CO2_moles * 3;
P1.O2_react_moles = P1.total_O_moles / 2;
P1.O2_mass = P1.O2_react_moles * C.M_O2; %kg

% Finding mass fractions
P1.react_mass_tot = P1.HDPE_mass + P1.O2_mass;
P1.HDPE_frac = P1.HDPE_mass / P1.react_mass_tot;
P1.O2_frac = P1.O2_mass / P1.react_mass_tot;

% we want to know the index of relevant species
P1.iH  = speciesIndex(P1.gas,'H');
P1.iH2  = speciesIndex(P1.gas,'H2');
P1.iO  = speciesIndex(P1.gas,'O');
P1.iO2  = speciesIndex(P1.gas,'O2');
P1.iOH  = speciesIndex(P1.gas,'OH');
P1.iC  = speciesIndex(P1.gas,'C');
P1.iCO  = speciesIndex(P1.gas,'CO');
P1.iCO2  = speciesIndex(P1.gas,'CO2');
P1.iH2O = speciesIndex(P1.gas,'H2O');
P1.iN2  = speciesIndex(P1.gas,'N2');
P1.iNO  = speciesIndex(P1.gas,'NO');
P1.iC2H4 = speciesIndex(P1.gas,'C2H4');

% Find adiabatic flame temperature
P1.mr = [linspace(1,1.75,20), linspace(1.75,3.5,50), linspace(3.5 ,10, 30)];
% mixture ratio = (mdot oxidizer/mdot fuel)

% mass fractions of reactants
P1.y_r = zeros(P1.nsp,length(P1.mr));
P1.y_r(P1.iC2H4,:) = 1./(1+P1.mr);
P1.y_r(P1.iO2,:) = P1.mr./(1+P1.mr);

P1.h_o_rxn = (P1.hForm_C2H4)./(1+P1.mr);
P1.h_o_HDPE = P1.h_o_rxn + P1.y_r(P1.iC2H4,:)*P1.dh_shift;
P1.hReactants = P1.HDPE_mass .* (P1.h_o_HDPE + P1.mr .* P1.hForm_O2);

for i = 1:length(P1.mr)
    set(P1.gas,'T',C.To,'P',C.P_chamber,'Y',P1.y_r(:,i));
    equilibrate(P1.gas,'HP');
    
    set(P1.gas,'P',C.P_chamber,'H',P1.hReactants(i));
    
    equilibrate(P1.gas,'HP');
    P1.T_adiabatic(i) = temperature(P1.gas);
    P1.rho(:,i) = density(P1.gas);
    P1.entropy(:,i) = entropy_mass(P1.gas);
    
    P1.y_prod(:,i) = massFractions(P1.gas);
    P1.X_prod(:,i) = moleFractions(P1.gas);
    P1.h_o(:,i) = enthalpy_mass(P1.gas);
end

%% 

% solving for gas information at the throat

P1.Mt_1 = 0.75; P1.Mt_2 = 0.85; P1.Mt_3 = 0.9;
for i = 1:length(P1.mr)
    P1.pres_throat_eq(i) = C.P_chamber;
    set(P1.gas, 'T', P1.T_adiabatic(i), 'P', P1.pres_throat_eq(i), 'Y', P1.y_prod(:,i));
    P1.Mt_eq = 0;
    
   while P1.Mt_eq < 1          %(P1.entropy(:,i) < P1.entropy_throat)
       if P1.Mt_eq < P1.Mt_1
           P1.dP = 1e4 ;
       elseif P1.Mt_eq >= P1.Mt_1 && P1.Mt_eq < P1.Mt_2
           P1.dP = 1e4;
       elseif P1.Mt_eq >= P1.Mt_2 && P1.Mt_eq < P1.Mt_3
           P1.dP = 1e2;
       else 
           P1.dP = 5e1;
       end 
      P1.pres_throat_eq(i) =  P1.pres_throat_eq(i) - P1.dP;
      setState_SP(P1.gas, [P1.entropy(i),P1.pres_throat_eq(i)]);
      equilibrate(P1.gas,'SP');
      P1.a = soundspeed(P1.gas);
      P1.Vt_eq(i) = sqrt(2*(P1.h_o(:,i) - enthalpy_mass(P1.gas)));
      P1.Mt_eq = P1.Vt_eq(i)/P1.a;
   end
   P1.T_throat_eq(i) = temperature(P1.gas);
   P1.rho_throat_eq(i) = density(P1.gas);
   P1.a_throat_eq(i) = soundspeed(P1.gas); 
   P1.y_throat_eq(:,i) = massFractions(P1.gas);
   P1.x_throat_eq(:,i) = moleFractions(P1.gas);
   P1.k_eq(i) = cp_mass(P1.gas)/cv_mass(P1.gas);
end

for i = 1:length(P1.mr)
    P1.pres_throat_fro(i) = C.P_chamber;
    set(P1.gas, 'T', P1.T_adiabatic(i), 'P', P1.pres_throat_fro(i), 'Y', P1.y_prod(:,i));
    P1.Mt_fro = 0;
   
   while P1.Mt_fro < 1          %(P1.entropy(:,i) < P1.entropy_throat)
       if P1.Mt_fro < P1.Mt_1
           P1.dP = 1e4;
       elseif P1.Mt_fro >= P1.Mt_1 && P1.Mt_fro < P1.Mt_2
           P1.dP = 1e3;
       elseif P1.Mt_fro >= P1.Mt_2 && P1.Mt_fro < P1.Mt_3
           P1.dP = 1e2;
       else 
           P1.dP = 5e1;
       end
      P1.pres_throat_fro(i) =  P1.pres_throat_fro(i) - P1.dP;
      setState_SP(P1.gas, [P1.entropy(i),P1.pres_throat_fro(i)]);
      P1.a = soundspeed(P1.gas);
      P1.Vt_fro(i) = sqrt(2*(P1.h_o(:,i) - enthalpy_mass(P1.gas)));
      P1.Mt_fro = P1.Vt_fro(i)/P1.a;
   end
   P1.Mt_fro_test(i) = P1.Mt_fro;
   P1.T_throat_fro(i) = temperature(P1.gas);
   P1.rho_throat_fro(i) = density(P1.gas);
   P1.a_throat_fro(i) = soundspeed(P1.gas); 
   P1.y_throat_fro(:,i) = massFractions(P1.gas);
   P1.k_fro(i) = cp_mass(P1.gas)/cv_mass(P1.gas);
end

P1.c_star_eq = C.P_chamber./(P1.rho_throat_eq .* P1.a_throat_eq);
P1.c_star_fro = C.P_chamber./(P1.rho_throat_fro .* P1.a_throat_fro);

%% PART 2 Return of the temperature calculations

P2.Pe = C.P_atm; % for a pefectly expanded nozzle  
P2.gas = P1.gas; 
P2.entropy = P1.entropy; % still isentropic yo

for i = 1:length(P1.mr)
    set(P2.gas, 'S', P2.entropy(i), 'P', P2.Pe, 'Y', P1.y_throat_eq(:,i));
    equilibrate(P2.gas,'SP');
    P2.Ve_eq(i) = sqrt(2*(P1.h_o(:,i) - enthalpy_mass(P2.gas)));
    P2.Te_eq(i) = temperature(P2.gas);
    P2.rho_exit_eq(i) = density(P2.gas);    
    P2.X_exit_eq(:,i) = moleFractions(P2.gas);
end
P2.epsilon_eq = (P1.rho_throat_eq .* P1.Vt_eq) ./ (P2.rho_exit_eq .* P2.Ve_eq);
P2.Cf_eq = P2.Ve_eq ./ P1.c_star_eq;

for i = 1:length(P1.mr)
    set(P2.gas, 'S', P2.entropy(i), 'P', P2.Pe, 'X', P1.X_prod(:,i));
    P2.Ve_fro(i) = sqrt(2*(P1.h_o(:,i) - enthalpy_mass(P2.gas)));
    P2.Te_fro(i) = temperature(P2.gas);
    P2.rho_exit_fro(i) = density(P2.gas);
    P2.X_exit_fro(:,i) = moleFractions(P2.gas);
end
P2.epsilon_fro = (P1.rho_throat_fro .* P1.Vt_fro) ./ (P2.rho_exit_fro .* P2.Ve_fro);
P2.Cf_fro = P2.Ve_fro ./ P1.c_star_fro;

%% Include experimental data

P3 = load('C:\Users\Ian\Box\ME140 Team Rocket\Project 7\Test analysis\testfire3.mat');
P3.A_throat = 1.824e-4; %m^2
P3.P_chamber = max(P3.chamP(P3.start_index:P3.final_index))*1e3; % [kPa] ---> [Pa]
P3.t1 = P3.time(P3.start_index);
P3.t2 = P3.time(P3.final_index);
P3.dt = P3.t2 - P3.t1;

P3.mfuel = P3.mfuel/1e3;
P3.mr = P3.m_O2 / P3.mfuel;

P3.mdot_fuel = P3.mfuel / P3.dt;
P3.mdot_O2 = P3.m_O2 / P3.dt;
P3.c_star_data = (P3.P_chamber * P3.A_throat) / (P3.mdot_fuel + P3.mdot_O2);

P3.mdot = P3.m_total / P3.dt;
P3.thrust = P3.thrust(P3.start_index:P3.final_index); 
P3.Ve_avg = mean(P3.thrust) / P3.mdot; % [m/s] thrust = ve  + pressure (assume pressure terms cancel)
P3.Cf = P3.Ve_avg / P3.c_star_data;

%
P4 = load('C:\Users\Ian\Box\ME140 Team Rocket\Project 7\Test analysis\testfire4.mat');
P4.A_throat = 1.824e-4; %m^2
P4.P_chamber = max(P4.chamP(P4.start_index:P4.final_index))*1e3; % [kPa] ---> [Pa]
P4.t1 = P4.time(P4.start_index);
P4.t2 = P4.time(P4.final_index);
P4.dt = P4.t2 - P4.t1;

P4.mfuel = P4.mfuel/1e3;
P4.mr = P4.m_O2 / P4.mfuel;

P4.mdot_fuel = P4.mfuel / P4.dt;
P4.mdot_O2 = P4.m_O2 / P4.dt;
P4.c_star_data = (P4.P_chamber * P4.A_throat) / (P4.mdot_fuel + P4.mdot_O2);

P4.mdot = P4.m_total / P4.dt;
P4.thrust = P4.thrust(P4.start_index:P4.final_index); 
P4.Ve_avg = mean(P4.thrust) / P4.mdot; % [m/s] thrust = ve  + pressure (assume pressure terms cancel)
P4.Cf = P4.Ve_avg / P4.c_star_data;

%
P5 = load('C:\Users\Ian\Box\ME140 Team Rocket\Project 7\Test analysis\testfire5.mat');
P5.A_throat = 1.824e-4; %m^2
P5.P_chamber = max(P5.chamP(P5.start_index:P5.final_index))*1e3; % [kPa] ---> [Pa]
P5.t1 = P5.time(P5.start_index);
P5.t2 = P5.time(P5.final_index);
P5.dt = P5.t2 - P5.t1;

P5.mfuel = P5.mfuel/1e3;
P5.mr = P5.m_O2 / P5.mfuel;

P5.mdot_fuel = P5.mfuel / P5.dt;
P5.mdot_O2 = P5.m_O2 / P5.dt;
P5.c_star_data = (P5.P_chamber * P5.A_throat) / (P5.mdot_fuel + P5.mdot_O2);

P5.mdot = P5.m_total / P5.dt;
P5.thrust = P5.thrust(P5.start_index:P5.final_index); 
P5.Ve_avg = mean(P5.thrust) / P5.mdot; % [m/s] thrust = ve  + pressure (assume pressure terms cancel)
P5.Cf = P5.Ve_avg / P5.c_star_data;

%
P6 = load('C:\Users\Ian\Box\ME140 Team Rocket\Project 7\Test analysis\testfire6.mat');
P6.A_throat = 1.824e-4; %m^2
P6.P_chamber = max(P6.chamP(P6.start_index:P6.final_index))*1e3; % [kPa] ---> [Pa]
P6.t1 = P6.time(P6.start_index);
P6.t2 = P6.time(P6.final_index);
P6.dt = P6.t2 - P6.t1;

P6.mfuel = P6.mfuel/1e3;
P6.mr = P6.m_O2 / P6.mfuel;

P6.mdot_fuel = P6.mfuel / P6.dt;
P6.mdot_O2 = P6.m_O2 / P6.dt;
P6.c_star_data = (P6.P_chamber * P6.A_throat) / (P6.mdot_fuel + P6.mdot_O2);

P6.mdot = P6.m_total / P6.dt;
P6.thrust = P6.thrust(P6.start_index:P6.final_index); 
P6.Ve_avg = mean(P6.thrust) / P6.mdot; % [m/s] thrust = ve  + pressure (assume pressure terms cancel)
P6.Cf = P6.Ve_avg / P6.c_star_data;

%
P7 = load('C:\Users\Ian\Box\ME140 Team Rocket\Project 7\Test analysis\testfire7.mat');
P7.A_throat = 1.824e-4; %m^2
P7.start_index = 2280; P7.final_index = 5250;
P7.P_chamber = max(P7.chamP(P7.start_index:P7.final_index))*1e3 + C.P_atm; % [kPa] ---> [Pa]
P7.t1 = P7.time(P7.start_index);
P7.t2 = P7.time(P7.final_index);
P7.dt = P7.t2 - P7.t1;

P7.mfuel = P7.mfuel/1e3;
P7.mr = P7.m_O2 / P7.mfuel;

P7.mdot_fuel = P7.mfuel / P7.dt;
P7.mdot_O2 = P7.m_O2 / P7.dt;
P7.c_star_data = (P7.P_chamber * P7.A_throat) / (P7.mdot_fuel + P7.mdot_O2);

P7.mdot = P7.m_total/P7.dt;
P7.thrust = P7.thrust(P7.start_index:P7.final_index); 
P7.Ve_avg = mean(P7.thrust)/P7.mdot; % [m/s] thrust = ve  + pressure (assume pressure terms cancel)
% V = M*c -> c = V/M = sqrt(krt)
P7.Cf = P7.Ve_avg / P7.c_star_data;

%% All the plots
% plot Cf vs mr
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
plot(P1.mr, P2.Cf_eq, '-', P1.mr, P2.Cf_fro, '-');
plot(P3.mr, P3.Cf, 'x');
plot(P4.mr, P4.Cf, 's');
plot(P5.mr, P5.Cf, '*');
plot(P6.mr, P6.Cf, 'd');
plot(P7.mr, P7.Cf, 'o');
legend('Equilibrium', 'Frozen', 'test 3', 'test 4', 'test 5', 'test 6', 'test 7', 'Location', 'best');
xlabel('Mixing Ratio, mr'); ylabel('C_F');
plotfixer;

% plot C* vs mr
figure('units','normalized'); % for larger plot
hold on;
plot(P1.mr, P1.c_star_eq, P1.mr, P1.c_star_fro);
plot(P3.mr, P3.c_star_data, 'x');
plot(P4.mr, P4.c_star_data, 's');
plot(P5.mr, P5.c_star_data, '*');
plot(P6.mr, P6.c_star_data, 'd');
plot(P7.mr, P7.c_star_data, 'o');
legend('Equilibrium', 'Frozen', 'test 3', 'test 4', 'test 5', 'test 6', 'test 7', 'Location', 'best');
xlabel('Mixing Ratio, mr'); ylabel('c*');
plotfixer;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on
plot(P1.mr, P2.Ve_eq, P1.mr, P2.Ve_fro);
plot(P3.mr, P3.Ve_avg, 'x');
plot(P4.mr, P4.Ve_avg, 's');
plot(P5.mr, P5.Ve_avg, '*');
plot(P6.mr, P6.Ve_avg, 'd');
plot(P7.mr, P7.Ve_avg, 'o');
legend('Equilibrium', 'Frozen', 'test 3', 'test 4', 'test 5', 'test 6', 'test 7', 'Location', 'best');
xlabel('Mixing Ratio, mr'); ylabel('Exit Velocity [m/s]');
plotfixer;

% figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
% plot(P1.mr, P2.Te_eq, P1.mr, P2.Te_fro);
% legend('Equilibrium','Frozen', 'Location', 'best');
% xlabel('Mixing Ratio, mr'); ylabel('T_{exit} [K]');
% plotfixer;

%% Export Files
E.mr = P1.mr;
E.mr3 = P3.mr; E.mr4 = P4.mr;
E.mr5 = P5.mr; E.mr6 = P6.mr;
E.mr7 = P7.mr;

E.Cf_eq = P2.Cf_eq; E.Cf_fro = P2.Cf_fro;
E.Cf3 = P3.Cf; E.Cf4 = P4.Cf;
E.Cf5 = P5.Cf; E.Cf6 = P6.Cf;
E.Cf7 = P7.Cf;
save('cFData', '-struct', 'E', 'mr', 'mr3', 'mr4', 'mr5', 'mr6', 'mr7', ...
    'Cf_eq', 'Cf_fro', 'Cf3', 'Cf4', 'Cf5', 'Cf7');

E.c_star_eq = P1.c_star_eq; E.c_star_fro = P1.c_star_fro;
E.c_star3 = P3.c_star_data; E.c_star4 = P4.c_star_data;
E.c_star5 = P5.c_star_data; E.c_star6 = P6.c_star_data;
E.c_star7 = P7.c_star_data;
save('cStarData', '-struct', 'E', 'mr', 'mr3', 'mr4', 'mr5', 'mr6', 'mr7', ...
    'c_star_eq', 'c_star_fro', 'c_star3', 'c_star4', 'c_star5',  'c_star6', 'c_star7');

E.Ve_eq = P2.Ve_eq; E.Ve_fro = P2.Ve_fro;
E.Ve3 = P3.Ve_avg; E.Ve4 = P4.Ve_avg;
E.Ve5 = P5.Ve_avg; E.Ve6 = P6.Ve_avg; 
E.Ve7 = P7.Ve_avg;
save('vExitData', '-struct', 'E', 'mr', 'mr3', 'mr4', 'mr5', 'mr6', 'mr7', ...
    'Ve_eq', 'Ve_fro', 'Ve3', 'Ve4', 'Ve5', 'Ve6', 'Ve7');

%% Impulse 
% assume correctly expanded and stuff so we get: F_T = m_dot * Ve
% P4.m_dot_eq = P2.rho_exit_eq .*P2.Ve_eq .* P3.A_throat .* P2.epsilon_eq;
% P4.F_t_eq = P4.m_dot_eq .* P2.Ve_eq;
% %P4.F_t_eq =  P4.m_dot_eq .* P1.c_star_eq .* P2.Cf_eq;
% 
% P4.m_dot_fro = P2.rho_exit_fro .* P2.Ve_fro .* P3.A_throat .* P2.epsilon_fro;
% P4.F_t_fro = P4.m_dot_fro .* P2.Ve_fro;
% %P4.F_t_fro = P4.m_dot_fro .*P1.c_star_fro .* P2.Cf_fro; 
% 
% figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
% plot(P1.mr, P4.m_dot_eq, P1.mr, P4.m_dot_fro)
% xlabel('Mixture Ratio, mr'); ylabel('mass flow');
% 
% figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
% plot(P1.mr, P4.F_t_eq, P1.mr, P4.F_t_fro);
% xlabel('Mixture Ratio, mr'); ylabel('Thrust, [N]');
% legend('Equilibrium','Frozen');
% plotfixer;

toc;
disp('finished running');
