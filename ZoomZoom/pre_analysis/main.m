%% Main Script
% ME 140 - Reginald Mitchell - Spring 2018
% Projects 7 - The Hybrid Rocket Motor
% Alex Casa, Ian Gunady, Minh Ngo, Thomas Trzpit, Zane Zook

% x is mol fraction
% y is mass fraction
%% Constants
clear all; close all; clc; 
tic

% C struct - constants
C.P_atm = 100000;           % [Pa] atmospheric pressure, 1 bar
C.P_chamber = 10*100000;    % [Pa] combustion chamber pressure, 68 bar
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
P1.HDPE_mass = 1; %kg
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

%% All the plots
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P1.c_star_eq = C.P_chamber./(P1.rho_throat_eq .* P1.a_throat_eq);
P1.c_star_fro = C.P_chamber./(P1.rho_throat_fro .* P1.a_throat_fro);
plot(P1.mr, P1.c_star_eq, P1.mr, P1.c_star_fro);
legend('Equilibrium', 'Frozen', 'Location', 'best');
xlabel('Mixing Ratio, mr'); ylabel('c*');
plotfixer;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P1.mr, P1.T_throat_eq, P1.mr, P1.T_throat_fro, P1.mr, P1.T_adiabatic);
xlabel('Mixing Ratio, mr'); ylabel('Temperature, [K]');
legend('Throat Equilibirum' , 'Throat Frozen' , 'Adiabatic Flame', 'Location', 'best');
plotfixer;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
plot(P1.mr,P1.X_prod(P1.iH,:));
plot(P1.mr,P1.X_prod(P1.iH2,:));
plot(P1.mr,P1.X_prod(P1.iO,:));
plot(P1.mr,P1.X_prod(P1.iO2,:));
plot(P1.mr,P1.X_prod(P1.iOH,:));
plot(P1.mr,P1.X_prod(P1.iC,:));
plot(P1.mr,P1.X_prod(P1.iCO,:));
plot(P1.mr,P1.X_prod(P1.iCO2,:));
plot(P1.mr,P1.X_prod(P1.iH2O,:));
plot(P1.mr,P1.X_prod(P1.iN2,:));
plot(P1.mr,P1.X_prod(P1.iNO,:));
plot(P1.mr,P1.X_prod(P1.iC2H4,:));
xlabel('Mixing Ratio, mr'); ylabel('Mole Fraction');
title('Frozen Mole Fractions at Combustor Exit')
legend('H', 'H_{2}' , 'O', 'O_{2}', 'OH', 'C', 'CO', 'CO_{2}', 'H_{2}O',...
    'N_2', 'NO', 'CH2', 'Location', 'best');
axis([1 10 0 1]);
plotfixer;


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


%% All the plots

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P1.mr, P2.Ve_eq, P1.mr, P2.Ve_fro);
legend('Equilibrium', 'Frozen', 'Location', 'best');
xlabel('Mixing Ratio, mr'); ylabel('Exit Velocity [m/s]');
plotfixer;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P1.mr, P2.Te_eq, P1.mr, P2.Te_fro);
legend('Equilibrium', 'Frozen', 'Location', 'best');
xlabel('Mixing Ratio, mr'); ylabel('T_{exit} [K]');
plotfixer;
%%
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P1.mr, P2.Cf_eq, '-', P1.mr, P2.Cf_fro, '-');
legend('Equilibrium', 'Frozen', 'Location', 'best');
xlabel('Mixing Ratio, mr'); ylabel('C_F');
plotfixer;
%%
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P1.mr, P2.epsilon_eq, P1.mr, P2.epsilon_fro);
legend('Equilibrium', 'Frozen', 'Location', 'best');
xlabel('Mixing Ratio, mr'); ylabel('Expansion Ratio');
plotfixer;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
plot(P1.mr, P2.X_exit_eq(P1.iH,:));
plot(P1.mr, P2.X_exit_eq(P1.iH2,:));
plot(P1.mr, P2.X_exit_eq(P1.iO,:));
plot(P1.mr, P2.X_exit_eq(P1.iO2,:));
plot(P1.mr, P2.X_exit_eq(P1.iOH,:));
plot(P1.mr, P2.X_exit_eq(P1.iC,:));
plot(P1.mr, P2.X_exit_eq(P1.iCO,:));
plot(P1.mr, P2.X_exit_eq(P1.iCO2,:));
plot(P1.mr, P2.X_exit_eq(P1.iH2O,:));
plot(P1.mr, P2.X_exit_eq(P1.iN2,:));
plot(P1.mr, P2.X_exit_eq(P1.iNO,:));
plot(P1.mr, P2.X_exit_eq(P1.iC2H4,:));
xlabel('Mixing Ratio, mr'); ylabel('Equilibrium Mole Fraction at Nozzle Exit');
legend('H', 'H_{2}' , 'O', 'O_{2}', 'OH', 'C', 'CO', 'CO_{2}', 'H_{2}O',...
    'N_2', 'NO', 'CH2', 'Location', 'EastOutside');
axis([1 10 0 1]);
plotfixer;

%%  P1.x_throat_eq
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
plot(P1.mr, P1.x_throat_eq(P1.iH,:));
plot(P1.mr, P1.x_throat_eq(P1.iH2,:));
plot(P1.mr, P1.x_throat_eq(P1.iO,:));
plot(P1.mr, P1.x_throat_eq(P1.iO2,:));
plot(P1.mr, P1.x_throat_eq(P1.iOH,:));
plot(P1.mr, P1.x_throat_eq(P1.iC,:));
plot(P1.mr, P1.x_throat_eq(P1.iCO,:));
plot(P1.mr, P1.x_throat_eq(P1.iCO2,:));
plot(P1.mr, P1.x_throat_eq(P1.iH2O,:));
plot(P1.mr, P1.x_throat_eq(P1.iN2,:));
plot(P1.mr, P1.x_throat_eq(P1.iNO,:));
plot(P1.mr, P1.x_throat_eq(P1.iC2H4,:));
xlabel('Mixing Ratio, mr'); ylabel('Throat Mole Fraction');
legend('H', 'H_{2}' , 'O', 'O_{2}', 'OH', 'C', 'CO', 'CO_{2}', 'H_{2}O',...
    'N_2', 'NO', 'CH2', 'Location', 'EastOutside');
axis([1 10 0 1]);
plotfixer;

%% Plot experimental data
P3 = load('RocketTestFire.mat');
P3.A_throat = 1.824e-4; %m^2
P3.P_chamber = mean(P3.chamP(P3.start_index:P3.final_index))*1e3; % [kPa] ---> [Pa]
P3.dt = (P3.final_index - P3.start_index)/1e3;

P3.m_fuel = P3.m_total - P3.m_O2;
P3.m_dot_fuel = P3.m_fuel/P3.dt;
P3.m_dot_O2 = mean(P3.m_dot_O2(P3.start_index:P3.final_index));

P3.mr = P3.m_dot_O2 / P3.m_dot_fuel;
P3.c_star_data = (P3.P_chamber * P3.A_throat) / (P3.m_dot_fuel + P3.m_dot_O2);

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
plot(P1.mr, P1.c_star_eq, P1.mr, P1.c_star_fro);
plot(P3.mr, P3.c_star_data, '*');
legend('Equilibrium', 'Frozen', 'Experimental Data', 'Location', 'best');
xlabel('Mixing Ratio, mr'); ylabel('c*');
plotfixer;

%% Impulse 
% assume correctly expanded and stuff so we get: F_T = m_dot * Ve
P4.m_dot_eq = P2.rho_exit_eq .*P2.Ve_eq .* P3.A_throat .* P2.epsilon_eq;
P4.F_t_eq = P4.m_dot_eq .* P2.Ve_eq;
%P4.F_t_eq =  P4.m_dot_eq .* P1.c_star_eq .* P2.Cf_eq;

P4.m_dot_fro = P2.rho_exit_fro .* P2.Ve_fro .* P3.A_throat .* P2.epsilon_fro;
P4.F_t_fro = P4.m_dot_fro .* P2.Ve_fro;
%P4.F_t_fro = P4.m_dot_fro .*P1.c_star_fro .* P2.Cf_fro; 

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P1.mr, P4.m_dot_eq, P1.mr, P4.m_dot_fro)
xlabel('Mixture Ratio, mr'); ylabel('mass flow');

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P1.mr, P4.F_t_eq, P1.mr, P4.F_t_fro);
xlabel('Mixture Ratio, mr'); ylabel('Thrust, [N]');
legend('Equilibrium','Frozen');
plotfixer;

toc;
disp('finished running');
