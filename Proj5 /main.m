%% Main Script
% ME 140 - Reginald Mitchell - Spring 2018
% Projects 5&6 - PEM Fuel Cell System Evaluation and Hydrogen Production Analysis
% Ian Gunady, Hayden Hall, Minh Ngo, Nick Wilson

%% Set up
close all; clear all; clc;
filename = 'PEMData.xlsx';
data = xlsread(filename);

% Constants
C.C2K = 273.15;                     % [K]
C.atm2Pa = 101400;                  % [Pa]
C.ft2m = 0.3048;                    % [m/ft]
C.R_u = 8.314;                      % [J/ mol K]
C.psi2Pa = 6894.76;
C.Patm = 1 *C.atm2Pa;
C.T_stand = 25 + C.C2K;
C.M_H2 = 2.01588e-3;                % [kg/mol]
C.M_H2O = 18.01528e-3;              % [kg/mol]
C.M_CO = 28.01e-3;                  % [kg/mol]
C.M_CO2 = 44.01e-3;                 % [kg/mol]
C.M_air = 28.9647e-3;               % [kg/mol]
C.stack_mult = 3;                   % some multiplier

% LHV, HHV of H2 from Table A-27
C.LHV_H2 = 120e6;               % [J/kg]
C.LHV_H2_bar = 120e6 * C.M_H2;  % [J/mol]

% Coefficients for specific heats
C.coef_H2Ov = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9];      %  specific heats in [J/mol]
C.coef_H2 = [29.11, -0.1916e-2, 0.4003e-5, -0.8704e-9];
C.coef_O2 = [25.48, 1.520e-2, -0.7155e-5, 1.312e-9];
C.coef_N2 = [28.90, -0.1571e-2, 0.8081e-5, -2.873e-9];
C.coef_CO = [28.16, 0.1675e-2, 0.5372e-5, -2.222e-9];
C.coef_CO2 = [22.26, 5.981e-2, -3.501e-5, 7.469e-9];
C.coef_CH4 = [19.89, 5.024e-2, 1.269e-5, -11.01e-9];

% Heat of formation and entropy of formation at STP

% N2 
C.h_fo_N2 = 0;                                              % [J/mol]
C.s_o_N2 = 191.61;                                          % [J/mol K]

% O2
C.h_fo_O2 = 0;                                              % [J/mol]
C.s_o_O2 = 205.04;                                          % [J/mol K] 

% H2
C.h_fo_H2 = 0;                                              % [J/mol]
C.s_o_H2 = 130.68;                                          % [J/mol K]

% H2O vapor 
C.h_fo_H2Ov = -241820;                                      % [J/mol]
C.s_o_H2Ov = 188.83;                                        % [J/mol K]

% H2O liquid
C.h_fo_H2Ol = -285830;                                      % [J/mol]
C.s_o_H2Ol = 69.92;                                         % [J/mol K]
C.cp_H2Ol = 4.18e3;                                         % [J/kg]
C.cp_bar_H2Ol = 4.18e3*C.M_H2O;                             % [J/mol K]

% CO (SMR, WGS)
C.h_fo_CO = -110530;                                        % [J/mol]
C.s_o_CO = 197.65;                                          % [J/mol K]

% CO2 (WGS)
C.h_fo_CO2 = -393520;                                       % [J/mol]
C.s_o_CO2 = 213.8;                                          % [J/mol K]

% CH4 (WGS)
C.h_fo_CH4 = -74850;                                        % [J/mol]
C.s_o_CH4 = 186.16;                                         % [J/mol K]

% Number of mols of fuel per reaction
C.N_fuel = 1;                                               % H2 

% Measured Data
D.NumR = data(1,:);
D.Resistor = data(2,:);
D.V_stack = data(3,:);                                      % [V]
D.I_stack = data(4,:) * C.stack_mult ;                      % [A]
D.V_load = data(5,:);                                       % [V]
D.I_load = data(6,:);                                       % [A]
D.V_dot_H2 =  data(7,:) * C.ft2m^3;                         % [SCFPH] --> [SCMPH]
D.V_dot_air = data(8,:) * C.ft2m^3;                         % [SCFPM] --> [SCMPM]
D.T1 = data(9,:) + C.C2K                                   % [K]
D.T2 = data(10,:) + C.C2K;
D.T3 = data(11,:) + C.C2K;
D.T4 = data(12,:) + C.C2K;
D.T5 = data(13,:) + C.C2K;
D.T6 = data(14,:) + C.C2K;
D.P_air = data(15,:) * C.psi2Pa +  C.Patm;  % [Pa]
D.P_H2 = data(16,:) *C.psi2Pa +  C.Patm;    % [Pa]

% Total resistance
D.R_tot(1:2) = D.Resistor(1:2);
for i = 3:length(D.Resistor)
    D.R_tot(i) = (1/D.R_tot(i-1) + 1/D.Resistor(i))^-1;
end 

%% Part A Problem 1 
% Power_stack = everything, all the power 
P1.P_stack = D.V_stack .* D.I_stack;                % [W]
% Power_load = "useful work", through the resistors 
P1.P_load = D.V_load .* D.I_load;                   % [W]
% Power_accessory = P_stack - P_load
P1.P_acc = P1.P_stack  - P1.P_load;                 % [W]

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P1.plot1 = plot(P1.P_load, D.V_stack, 'o-', P1.P_load, D.V_load, 'o-');
xlabel('Load Power'); ylabel('Voltage, [V]');
plotfixer;
legend('Stack Voltage', 'Load Voltage', 'Location', 'best');

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P1.plot2 = plot(P1.P_load, D.I_stack, 'o-', P1.P_load, D.I_load, 'o-');
xlabel('Load Power'); ylabel('Current, [A]');
plotfixer;
legend('Stack Current', 'Load Current', 'Location', 'best');

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P1.plot3 = plot(P1.P_load, P1.P_stack, 'o-', P1.P_load, P1.P_acc, 'o-');
xlabel('Load Power'); ylabel('Power, [W]');
plotfixer;
legend('Stack Power', 'Accessory Power', 'Location', 'best');

% [SCMPM] -> [mol/sec] convert Standard Cubic Meters Per Minute to mol/sec
P1.Psat_std = T2P_sat(C.T_stand);
P1.P_dryAir = C.Patm - 0.36 * P1.Psat_std;
P1.n_dot_air = P1.P_dryAir * D.V_dot_air /(C.R_u * C.T_stand)/60;
P1.m_dot_air = P1.n_dot_air *C.M_air; 

% [SCMFH] -> [mol/sec]  convert Standard Cubic Meters Per Hour to mol/sec
P1.n_dot_H2 = C.Patm * D.V_dot_H2 /(C.R_u * C.T_stand)/3600;
P1.m_dot_H2 = P1.n_dot_H2 *C.M_H2; 

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P1.plot3 = plot(P1.P_load, P1.m_dot_air*1e3 , 'o-', P1.P_load, P1.m_dot_H2*1e5, 'o-');
xlabel('Load Power, [W]'); ylabel('Mass Flow Rate,[g/s]');
plotfixer;
legend('Air Flow', '100 x H_2 Flow', 'Location', 'best')

%% Part A, Problem 2
% Reactants: H2 + 0.5*L*(O2+3.76*N2) + alpha*H2Ov
% Products: beta*H2Ov + gamma*H2Ol + 0.5*(L-1)*O2 + 0.5*L*3.76*N2

P2.Pm = C.Patm;                             % [Pa] 
P2.n_dot_O2 = 1/(1+3.76) * P1.n_dot_air;    % [mol/sec]
P2.L = 2*P2.n_dot_O2 ./ P1.n_dot_H2;        % [unitless]

P2.Psat_react = T2P_sat(D.T1);              % [Pa]
P2.Psat_prod = T2P_sat(D.T2);               % [Pa]

P2.y_max_react = P2.Psat_react./P2.Pm;      % [unitless]
P2.y_max_prod = P2.Psat_prod./P2.Pm;        % [unitless]

% 100% RH, saturation at the inlet 
P2.alpha = P2.y_max_react./(1-P2.y_max_react).*(0.5*P2.L*(1+ 3.76));

P2.N_H2O_prod = 1 + P2.alpha;
P2.N_a_prod =  0.5*(P2.L - 1) + 0.5*P2.L*(3.76);

P2.y_test_prod = P2.N_H2O_prod ./ (P2.N_H2O_prod + 0.5*(P2.L-1) + 0.5*P2.L*3.76); 

[P2.beta, P2.gamma] = vaporLiquidBalance (D.T2, P2.y_test_prod, ...
                        P2.y_max_prod, P2.N_a_prod, P2.N_H2O_prod, P2.alpha);

P2.N_prod = P2.beta + P2.N_a_prod;
P2.N_react = 0.5*(P2.L)*(1+3.76) + P2.alpha;

% N2 REACTANT
P2.N_N2_react = (.5*P2.L*3.76);                                     
P2.y_N2_react = P2.N_N2_react./(P2.N_react);                        % mol fraction of N2 on reactants side
P2.prat_N2_react = P2.y_N2_react.*D.P_air./C.Patm;                  
% pressure ratio
P2.ds_o_N2_react = delta_s (C.coef_N2, D.T1, P2.prat_N2_react);     
% integral term and log term of entropy 
P2.dh_bar_N2_react = integral_h (C.coef_N2, D.T1);                  % [J/mol]

P2.g_bar_N2react = C.h_fo_N2 + P2.dh_bar_N2_react ...
                    - D.T1.*(C.s_o_N2 + P2.ds_o_N2_react);          % [J/mol]

% N2 PRODUCT 
P2.N_N2_prod = (.5*P2.L*3.76);
P2.y_N2_prod = P2.N_N2_prod./P2.N_prod;                             
% mole fraction of N2 on PRODUCT side
P2.ds_o_N2_prod = delta_s (C.coef_N2, D.T2, P2.y_N2_prod);          
% integral term and log term of entropy 
P2.dh_bar_N2_prod = integral_h (C.coef_N2, D.T2);                   % [J/mol]

P2.g_bar_N2_prod = C.h_fo_N2 + P2.dh_bar_N2_prod ...
                    - D.T2.*(C.s_o_N2 + P2.ds_o_N2_prod);

% O2 REACTANT 
P2.N_O2_react = (.5*P2.L); 
P2.y_O2_react = P2.N_O2_react ./(P2.N_react);                       
% mole fraction of O2 on reactants side
P2.prat_O2_react = P2.y_O2_react.*D.P_air./C.Patm;                  
% pressure ratio
P2.ds_o_O2_react = delta_s (C.coef_O2, D.T1, P2.prat_O2_react);     
% integral term and log term of entropy 
P2.dh_bar_O2_react = integral_h (C.coef_O2, D.T1);       % [J/mol]

P2.g_bar_O2react = C.h_fo_O2 + P2.dh_bar_O2_react ...
                    - D.T1.*(C.s_o_O2 + P2.ds_o_O2_react);
        
% O2 PRODUCT 
P2.N_O2_prod = (0.5*(P2.L - 1));
P2.y_O2_prod = P2.N_O2_prod ./P2.N_prod;                            
% mole fraction of O2 on reactants side
P2.ds_o_O2_prod = delta_s (C.coef_O2, D.T2, P2.y_O2_prod);          
% integral term and log term of entropy 
P2.dh_bar_O2_prod = integral_h (C.coef_O2, D.T2);       % [J/mol]

P2.g_bar_O2_prod = C.h_fo_O2 + P2.dh_bar_O2_prod ...
                    - D.T2.*(C.s_o_O2 + P2.ds_o_O2_prod);

% H2 REACTANT
P2.prat_H2_react = D.P_H2./C.Patm;
P2.ds_o_H2 = delta_s (C.coef_H2, D.T1, P2.prat_H2_react);      
% integral term and log term of entropy
P2.dh_bar_H2 = integral_h (C.coef_H2, D.T1);        % [J/mol]

P2.g_bar_H2 = C.h_fo_H2 + P2.dh_bar_H2 ...
                - D.T1.*(C.s_o_H2 + P2.ds_o_H2);

% H2O vapor REACTANT
P2.y_H2Ov_react = P2.alpha./P2.N_react;
P2.prat_H2Ov_react = P2.y_H2Ov_react.*D.P_air./C.Patm;
P2.ds_o_H2Ov_react = delta_s (C.coef_H2Ov, D.T1, P2.prat_H2Ov_react);        
% integral term and log term of entropy
P2.dh_bar_H2Ov_react = integral_h (C.coef_H2Ov, D.T1);    % [J/mol]

P2.g_bar_H2Ov_react = C.h_fo_H2Ov + P2.dh_bar_H2Ov_react ...
                        - D.T1.*(C.s_o_H2Ov + P2.ds_o_H2Ov_react);
        
% H2O vapor PRODUCT   
P2.y_H2Ov_prod = P2.beta./P2.N_prod; 
P2.ds_o_H2Ov_prod = delta_s (C.coef_H2Ov, D.T2, P2.y_H2Ov_prod);        
% integral term and log term of entropy
P2.dh_bar_H2Ov_prod = integral_h (C.coef_H2Ov, D.T2);    % [J/mol]

P2.g_bar_H2Ov_prod = C.h_fo_H2Ov + P2.dh_bar_H2Ov_prod ...
                    - D.T2.*(C.s_o_H2Ov + P2.ds_o_H2Ov_prod);

% H2O Liquid PRODUCT
P2.dh_bar_H2Ol_prod = C.cp_bar_H2Ol.*(D.T2 - C.T_stand);    % [J/mol]
P2.ds_o_H2Ol_prod = C.cp_bar_H2Ol.*log(D.T2/C.T_stand); 

P2.g_bar_H2Ol_prod = C.h_fo_H2Ol + P2.dh_bar_H2Ol_prod ...
                         - D.T2.*(C.s_o_H2Ol + P2.ds_o_H2Ol_prod);

% Delta G [J/mol H2]
P2.DG = P2.beta .* P2.g_bar_H2Ov_prod ...
        + P2.gamma .* P2.g_bar_H2Ol_prod ...
        + 0.5 *(P2.L - 1).* P2.g_bar_O2_prod ...
        + (.5*P2.L*3.76).* P2.g_bar_N2_prod...
        - P2.g_bar_H2...                                    % N_H2 = 1  
        - 0.5*P2.L.*P2.g_bar_O2react...                     % N_O2_react = 1 = 1/2 * 2 
        - (.5*P2.L.*3.76).* P2.g_bar_N2react...
        - P2.alpha.*P2.g_bar_H2Ov_react;     

P2.P_H2O_v = P2.alpha./(P2.alpha+.5.*P2.L.*(1+3.76)).*P2.Pm;
P2.RH_inlet = P2.y_H2Ov_react.*P2.Pm./P2.Psat_react;

P2.LHV = C.LHV_H2_bar * C.N_fuel;                            % [J]

% First Law Efficiency 
P2.E_I_load = P1.P_load./(C.LHV_H2.*P1.m_dot_H2);
P2.E_I_stack = P1.P_stack./(C.LHV_H2.*P1.m_dot_H2);

% Second Law Efficiency
P2.E_II_load = P1.P_load./(-P2.DG .* P1.n_dot_H2);
P2.E_II_stack = P1.P_stack./(-P2.DG .* P1.n_dot_H2);

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P2.plot1 = plot(P1.P_load, P2.L , 'o-');
xlabel('Load Power, [W]'); ylabel('Excess-air Coefficient \lambda');
plotfixer;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P2.plot2 = plot(P1.P_load, P2.E_I_load , 'o-', P1.P_load, P2.E_I_stack , 'o-');
hold on; 
P2.plot2 = plot(P1.P_load, P2.E_II_load , 'o-', ...
                P1.P_load, P2.E_II_stack , 'o-');
legend('Load 1st Efficiency','Stack 1st Efficiency ',...
    'Load 2nd Efficiency','Stack 2nd Efficiency ',...
    'Location', 'best');
xlabel('Load Power, [W]'); ylabel('First and Second Law Efficiency');
plotfixer;

% Calculate Lost Power in Stack Only and Full System
P2.power_loss_load = -P2.DG.* P1.n_dot_H2 - P1.P_load;
P2.power_loss_stack = -P2.DG.* P1.n_dot_H2 - P1.P_stack;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P2.plot4 = plot(P1.P_load, P2.power_loss_load, 'o-',...
                P1.P_load, P2.power_loss_stack, '*-');
xlabel('Load Power, [W]'); ylabel('Power Loss, [W]');
legend('Full System','Stack Only', 'Location', 'best')
plotfixer;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P2.plot4 = plot(P1.P_stack, P2.power_loss_load, 'o-', ...
                P1.P_stack, P2.power_loss_stack, '*-');
xlabel('Stack Power, [W]'); ylabel('Power Loss, [W]');
legend('Full System','Stack Only', 'Location', 'best')
plotfixer;

%% Part A, Problem 3
C.cp_H2Ol = 4.18e3;                             % [J/kg]
C.gpm2Ls = 0.06309;                             % [gpm/(L/s)]
P3.m_dot_cool = 0.35 * C.gpm2Ls * 1;            % [gpm] --> [kg/s]
P3.Q_dot_cool1 = P3.m_dot_cool * C.cp_H2Ol *(D.T1(5) - D.T2(5));
P3.P_peak = 597;                                %[W]
P3.P_net = 90e3;                                %[W]
P3.scale = P3.P_net/P3.P_peak;
P3.Q_dot_cool_scale = P3.scale * P3.Q_dot_cool1;
P3.m_dot_cool_scale = P3.m_dot_cool * P3.scale;

%% Part B, Problem 1
% Kp for SMR and WGS reaction
P4.T = linspace(25, 1200, 40) + C.C2K; 
P4.Kp_SMR = Kp_SMR_function(P4.T); 
P4.Kp_WGS = Kp_WGS_function(P4.T); 

%% A semi-log plot (i.e., plot log10 KP vs. T).
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
P4.plot1 = semilogy(P4.T, P4.Kp_SMR , '-', ...
                        P4.T, P4.Kp_WGS, '-');
hold on;
P4.plot1 = semilogy(P4.T, ones(length(P4.T)) ,'b:');
xlabel('Temperature, [K]'); ylabel('K_p');
axis([P4.T(1), P4.T(end), 1e-3, 1e3]);
legend('SMR','WGS', 'Location', 'best')
plotfixer;

%% Part B, Problem 2
% plot mol fractions of SMR reaction
P5.Prat = [1,10,100];
P5.T = P4.T;

% CH4 + 3*H2O -> x*CO + y*H2 + z*CH4 + w*H2O 
% CH4 + H2O -> CO +3yH2
P5.Kp = P4.Kp_SMR;
P5.dx = 1e-5; 
P5.x = zeros(length(P5.Prat), length(P5.Kp)); 
P5.Kp_est = zeros(length(P5.Prat), length(P5.Kp)); 

for j = 1:length(P5.Prat)
    for i = 1:length(P5.Kp)
        while P5.Kp(i) > P5.Kp_est(j,i)
            P5.y(j,i) = 3 * P5.x(j,i);          % N_H2
            P5.z(j,i) = 1 - P5.x(j,i);          % N_CH4
            P5.w(j,i) = 3 - P5.x(j,i);          % N_H2O
            
            if (P5.x(j,i) >= 0 && P5.y(j,i) >= 0 && P5.w(j,i)>=0 && P5.z(j,i) >= 0)
                P5.Kp_est(j,i) = (P5.x(j,i) * P5.y(j,i)^3)/...
                            (P5.w(j,i) * P5.z(j,i)) ...
                    * (P5.Prat(j)/(P5.x(j,i)+ P5.y(j,i) + P5.z(j,i) + P5.w(j,i)))^2;
            else 
                break
            end
            P5.x(j,i) = P5.x(j,i) + P5.dx; 
        end
    end 
end
P5.y_total = P5.x + P5.y + P5.z + P5.w;
P5.y_CO = P5.x ./P5.y_total;
P5.y_H2 = P5.y ./P5.y_total;
P5.y_CH4 = P5.z ./P5.y_total;
P5.y_H2O = P5.w ./P5.y_total;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P5.T, P5.y_CO(1,:),'-b' ,...
        P5.T, P5.y_CO(2,:),'--b',...
        P5.T, P5.y_CO(3,:),'-.b')
hold on; 
plot( P5.T, P5.y_H2(1,:),'-r',...
        P5.T, P5.y_H2(2,:),'--r',...
        P5.T, P5.y_H2(3,:),'-.r')
    
plot( P5.T, P5.y_CH4(1,:),'-k',...
        P5.T, P5.y_CH4(2,:),'--k',...
        P5.T, P5.y_CH4(3,:),'-.k')
    
plot( P5.T, P5.y_H2O(1,:),'-m',...
        P5.T, P5.y_H2O(2,:),'--m',...
        P5.T, P5.y_H2O(3,:),'-.m');
    
axis([P5.T(1), P5.T(end), 0, 1]);
legend('CO 1atm', 'CO 10atm', 'CO 100atm',...
    'H2 1atm', 'H2 10atm', 'H2 100atm',...
    'CH4 1atm', 'CH4 10atm', 'CH4 100atm',...
    'H2O 1atm', 'H2O 10atm', 'H2O 100atm',...
    'Location', 'eastoutside');
ylabel('Mole Fraction, y'); xlabel('Temperature, [K]');
plotfixer;

%% Part B Problem 3
% Stoichiometric: H2O + CO -> CO2 + H2
% WGS Reaction: H2O + CO -> xCO2 + yH2 + zH2O + wCO;
P6.T = P4.T;
P6.Kp = P4.Kp_WGS;
P6.dx = 1e-6;
P6.x = zeros(1,length(P6.Kp)); 
P6.Kp_est = zeros(1, length(P6.Kp)); 

for i = 1:length(P6.Kp)
    while P6.Kp(i) > P6.Kp_est(i)
        P6.y(i) = P6.x(i) + 3;          % N_H2
        P6.z(i) = 2 - P6.x(i);          % N_CH4
        P6.w(i) = 1 - P6.x(i);          % N_H2O

        if (P6.x(i) >= 0 && P6.y(i) >= 0 && P6.w(i)>=0 && P6.z(i) >= 0)
            P6.Kp_est(i) = (P6.x(i) * P6.y(i))/(P6.w(i) * P6.z(i));
        else 
            break
        end
        P6.x(i) = P6.x(i) + P6.dx; 
    end
end 

P6.y_total = 6;
P6.y_CO2 = P6.x ./ P6.y_total;
P6.y_H2 = P6.y ./ P6.y_total;
P6.y_CO = P6.w ./ P6.y_total;
P6.y_H2O = P6.z ./ P6.y_total;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P6.T, P6.y_CO,'-b') 
hold on; 
plot( P6.T, P6.y_H2,'-r')
        
plot( P6.T, P6.y_CO2,'-k')
    
plot( P6.T, P6.y_H2O,'-m')
        
axis([P6.T(1), P6.T(end), 0, 0.7]);
legend('CO','H2','CO_2','H_2O',...
    'Location', 'best')
ylabel('Mole Fraction, y'); xlabel('Temperature, [K]');
plotfixer;

%% Part B Problem 4
% Station a
P7a.T_reform = 800 + C.C2K; 
P7a.Kp_reform = Kp_WGS_function(P7a.T_reform);

% WGS Reaction: H2O + CO -> xCO2 + yH2 + zH2O + wCO;

% Find composition of mixture from Kp; using fzero 
P7a.func = @(x) x*(x+3)/((1-x)*(2-x)) - P7a.Kp_reform;
P7a.x = fzero(P7a.func, 0);
P7a.y = P7a.x + 3;
P7a.z = 2 - P7a.x;
P7a.w = 1 - P7a.x;

% Find mol fractions
P7a.y_total = P7a.x + P7a.y + P7a.z + P7a.w;
P7a.y_CO2 = P7a.x ./P7a.y_total;
P7a.y_H2 = P7a.y ./P7a.y_total;
P7a.y_CO = P7a.w ./P7a.y_total;
P7a.y_H2O = P7a.z ./P7a.y_total;

% Change in enthalpy from STP
P7a.dh_bar_H2Ov = integral_h (C.coef_H2Ov, P7a.T_reform);       % [J/mol]
P7a.dh_bar_CH4 = integral_h (C.coef_CH4, P7a.T_reform);         % [J/mol]
P7a.dh_bar_CO2 = integral_h (C.coef_CO2, P7a.T_reform);         % [J/mol]
P7a.dh_bar_H2 = integral_h (C.coef_H2, P7a.T_reform);           % [J/mol]
P7a.dh_bar_CO = integral_h (C.coef_CO, P7a.T_reform);           % [J/mol]

% Consider everything coming in and then everything coming out.
P7a.DH_mol = P7a.y*(C.h_fo_H2 + P7a.dh_bar_H2) ...
        + P7a.x*(C.h_fo_CO2 + P7a.dh_bar_CO2)...
        + P7a.w*(C.h_fo_CO + P7a.dh_bar_CO)...
        + (P7a.z-3)*(C.h_fo_H2Ov + P7a.dh_bar_H2Ov)...
        - 1*(C.h_fo_CH4 + P7a.dh_bar_CH4);

P7a.m_react = (16.042 + 3*18.016)*1e-3; 
P7a.DH_kg = P7a.DH_mol /P7a.m_react /1e6;

C.MCH4 = 16.04e-3;      %[kg/mol]
C.LHV_CH4 = 50.05e6;        %[J/kg]
C.LHV_bar_CH4 = C.LHV_CH4 * C.MCH4;
P7a.burnedCH4_pc = P7a.DH_mol/(P7a.DH_mol + C.LHV_bar_CH4);

%% Station b
P7b.T_shift = 400 + C.C2K;
P7b.Kp_shift = Kp_WGS_function(P7b.T_shift);

% WGS Reaction: H2O + CO -> xCO2 + yH2 + zH2O + wCO;

% Find composition of mixture from Kp
P7b.func = @(x) x*(x+3)/((1-x)*(2-x)) - P7b.Kp_shift;
P7b.x = fzero(P7b.func, 0.5);
P7b.y = P7b.x + 3;
P7b.z = 2 - P7b.x;
P7b.w = 1 - P7b.x;

% Find mol fractions
P7b.y_total = P7b.x + P7b.y + P7b.z + P7b.w;
P7b.y_CO2 = P7b.x ./ P7b.y_total;
P7b.y_H2 = P7b.y ./ P7b.y_total;
P7b.y_CO = P7b.w ./ P7b.y_total;
P7b.y_H2O = P7b.z ./ P7b.y_total;

% Change in enthalpy from STP
P7b.dh_bar_H2Ov = integral_h (C.coef_H2Ov, P7b.T_shift);% [J/mol]
P7b.dh_bar_CO2 = integral_h (C.coef_CO2, P7b.T_shift);  % [J/mol]
P7b.dh_bar_H2 = integral_h (C.coef_H2, P7b.T_shift);    % [J/mol]
P7b.dh_bar_CO = integral_h (C.coef_CO, P7b.T_shift);    % [J/mol]

% Consider everything coming in and then everything coming out
P7b.DH_mol = (P7b.z-P7a.z)*(C.h_fo_H2Ov + P7b.dh_bar_H2Ov) ...
            + (P7b.w-P7a.w)*(C.h_fo_CO + P7b.dh_bar_CO)...
            + (P7b.x-P7a.x)*(C.h_fo_CO2 + P7b.dh_bar_CO2)...
            + (P7b.y-P7a.y)*(C.h_fo_H2 + P7b.dh_bar_H2);
        
P7b.m_react = P7a.z*C.M_H2O + P7a.w*C.M_CO + P7a.x*C.M_CO2 + P7a.y*C.M_H2;  
% [kg/mol(reaction)]
P7b.DH_kg = P7b.DH_mol/P7b.m_react/1e6;                 % [MJ/kg(reactant)]

%% Station c 

P7c.dT = 1e-1; 
P7c.T = P7b.T_shift;
P7c.DH_mol = P7b.DH_mol;
P7c.H_mol_react = (P7a.z)*(C.h_fo_H2Ov + P7b.dh_bar_H2Ov) ...
            + (P7a.w)*(C.h_fo_CO + P7b.dh_bar_CO)...
            + (P7a.x)*(C.h_fo_CO2 + P7b.dh_bar_CO2)...
            + (P7a.y)*(C.h_fo_H2 + P7b.dh_bar_H2);

while P7c.DH_mol < 0  
    P7c.T = P7c.T + P7c.dT; 
    P7c.Kp = Kp_WGS_function(P7c.T);
    P7c.func = @(x) x*(x+3)/((1-x)*(2-x)) - P7c.Kp;
    P7c.x = fzero(P7c.func, 0.5);
    
    % WGS Reaction: H2O + CO -> xCO2 + yH2 + zH2O + wCO;
    P7c.y = P7c.x + 3;
    P7c.z = 2 - P7c.x;
    P7c.w = 1 - P7c.x;

    P7c.y_total = 6;
    P7c.y_CO2 = P7c.x ./ P7c.y_total;
    P7c.y_H2 = P7c.y ./ P7c.y_total;
    P7c.y_CO = P7c.w ./ P7c.y_total;
    P7c.y_H2O = P7c.z ./ P7c.y_total;

    P7c.dh_bar_H2Ov = integral_h (C.coef_H2Ov, P7c.T);  % [J/mol]
    P7c.dh_bar_CO2 = integral_h (C.coef_CO2, P7c.T);    % [J/mol]
    P7c.dh_bar_H2 = integral_h (C.coef_H2, P7c.T);      % [J/mol]
    P7c.dh_bar_CO = integral_h (C.coef_CO, P7c.T);      % [J/mol]

    % Enthalpy at outlet - enthalpy at inlet
    P7c.DH_mol = (P7c.z)*(C.h_fo_H2Ov + P7c.dh_bar_H2Ov) ...
                + (P7c.w)*(C.h_fo_CO + P7c.dh_bar_CO)...
                + (P7c.x)*(C.h_fo_CO2 + P7c.dh_bar_CO2)...
                + (P7c.y)*(C.h_fo_H2 + P7c.dh_bar_H2)...
                - P7c.H_mol_react;
end

%% Station d
P7d.T_shift = 250 + C.C2K;
P7d.Kp_shift = Kp_WGS_function(P7d.T_shift);

P7d.dx = 1e-5; 
P7d.x = 0; 
P7d.Kp_est = 0; 

% WGS Reaction: H2O + CO -> x CO2 + y H2 +  z H2O + w CO;
while P7d.Kp_shift > P7d.Kp_est
    P7d.y = P7d.x +3;
    P7d.z = 2 - P7d.x;
    P7d.w = 1 - P7d.x;

    if (P7d.x >= 0 && P7d.y >= 0 && P7d.w>=0 && P7d.z >= 0)
        P7d.Kp_est = (P7d.x * P7d.y)/(P7d.w * P7d.z);
    else 
        break
    end
    P7d.x = P7d.x + P7d.dx; 
end

P7d.y_total = P7d.x + P7d.y + P7d.z + P7d.w;
P7d.y_CO2 = P7d.x ./P7d.y_total;
P7d.y_H2 = P7d.y ./P7d.y_total;
P7d.y_CO = P7d.w ./P7d.y_total;
P7d.y_H2O = P7d.z ./P7d.y_total;

P7d.dh_bar_H2Ov = integral_h (C.coef_H2Ov, P7d.T_shift);    % [J/mol]
P7d.dh_bar_CO2 = integral_h (C.coef_CO2, P7d.T_shift);      % [J/mol]
P7d.dh_bar_H2 = integral_h (C.coef_H2, P7d.T_shift);        % [J/mol]
P7d.dh_bar_CO = integral_h (C.coef_CO, P7d.T_shift);        % [J/mol]

% Consider everything coming in and then everything coming out
P7d.DH_mol = (P7d.z-P7b.z)*(C.h_fo_H2Ov + P7d.dh_bar_H2Ov) ...
            + (P7d.w-P7b.w)*(C.h_fo_CO + P7d.dh_bar_CO)...
            + (P7d.x-P7b.x)*(C.h_fo_CO2 + P7d.dh_bar_CO2)...
            + (P7d.y-P7b.y)*(C.h_fo_H2 + P7d.dh_bar_H2);
        
P7d.m_react = P7b.z*C.M_H2O + P7b.w*C.M_CO + P7b.x*C.M_CO2 + P7b.y*C.M_H2;  % [kg/mol(reaction)]
P7d.DH_kg = P7d.DH_mol/P7d.m_react/1e6;                 % [MJ/kg(reactant)]

P7d.N_H2 = P7d.y;
P7d.Eff = P7d.N_H2 * C.LHV_H2_bar /(1*C.LHV_bar_CH4 + P7a.DH_mol);

%% Station e

P7e.dT = [1, 1e-1, 1e-2, 1e-3];
P7e.T = P7d.T_shift;
P7e.DH_mol = P7d.DH_mol;
P7e.DH_mol_react = (P7c.z)*(C.h_fo_H2Ov + P7d.dh_bar_H2Ov) ...
            + (P7c.w)*(C.h_fo_CO + P7d.dh_bar_CO)...
            + (P7c.x)*(C.h_fo_CO2 + P7d.dh_bar_CO2)...
            + (P7c.y)*(C.h_fo_H2 + P7d.dh_bar_H2);

while (P7e.DH_mol < 0 && P7e.T < 750)
    P7e.T = P7e.T + P7e.dT(2); 
    P7e.Kp = Kp_WGS_function(P7e.T);
    P7e.x = find_x_from_KpWGS(P7e.Kp);

    % WGS Reaction: H2O + CO -> xCO2 + yH2 + zH2O + wCO;

    P7e.y = P7e.x + 3;
    P7e.z = 2 - P7e.x;
    P7e.w = 1 - P7e.x;

    P7e.y_total = 6;
    P7e.y_CO2 = P7e.x ./ P7e.y_total;
    P7e.y_H2 = P7e.y ./ P7e.y_total;
    P7e.y_CO = P7e.w ./ P7e.y_total;
    P7e.y_H2O = P7e.z ./ P7e.y_total;

    P7e.dh_bar_H2Ov = integral_h (C.coef_H2Ov, P7e.T);      % [J/mol]
    P7e.dh_bar_CO2 = integral_h (C.coef_CO2, P7e.T);        % [J/mol]
    P7e.dh_bar_H2 = integral_h (C.coef_H2, P7e.T);          % [J/mol]
    P7e.dh_bar_CO = integral_h (C.coef_CO, P7e.T);          % [J/mol]

    % Consider everything coming in and then everything coming out
    P7e.DH_mol = (P7e.z)*(C.h_fo_H2Ov + P7e.dh_bar_H2Ov) ...
                + (P7e.w)*(C.h_fo_CO + P7e.dh_bar_CO)...
                + (P7e.x)*(C.h_fo_CO2 + P7e.dh_bar_CO2)...
                + (P7e.y)*(C.h_fo_H2 + P7e.dh_bar_H2)...
                - P7e.DH_mol_react;
end

P7e.N_H2 = P7e.y;
P7e.Eff = P7e.N_H2 * C.LHV_H2_bar /(1*C.LHV_bar_CH4 + P7a.DH_mol);

%% plotting section

P7.stations = [100, 200, 300];
P7.y_CO_isotherm = [P7a.y_CO, P7b.y_CO, P7d.y_CO];
P7.y_CO2_isotherm = [P7a.y_CO2, P7b.y_CO2, P7d.y_CO2];
P7.y_H2O_isotherm = [P7a.y_H2O, P7b.y_H2O, P7d.y_H2O];
P7.y_H2_isotherm = [P7a.y_H2, P7b.y_H2, P7d.y_H2];

P7.y_CO_adiab = [P7a.y_CO, P7c.y_CO, P7e.y_CO];
P7.y_CO2_adiab = [P7a.y_CO2, P7c.y_CO2, P7e.y_CO2];
P7.y_H2O_adiab = [P7a.y_H2O, P7c.y_H2O, P7e.y_H2O];
P7.y_H2_adiab = [P7a.y_H2, P7c.y_H2, P7e.y_H2];

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot

plot(P7.stations, P7.y_CO_isotherm, 'b-o', ...
        P7.stations, P7.y_CO_adiab, 'b:o',...
        P7.stations, P7.y_CO2_isotherm,'r-o', ...
        P7.stations, P7.y_CO2_adiab, 'r:o', ...
        P7.stations, P7.y_H2_isotherm, 'k-o', ...
        P7.stations, P7.y_H2_adiab, 'k:o' ,...
        P7.stations, P7.y_H2O_isotherm,'m-o', ...
        P7.stations, P7.y_H2O_adiab, 'm:o');
legend('CO isothermal','CO adiabatic','CO_2 isothermal','CO_2 adiabatic',...
    'H_2 isothermal','H_2 adiabatic','H_2O isothermal','H2_O adiabatic')
set(gca,'XTickLabel',{'Reformer', 'Reactor 1', 'Reactor 2'})
xlabel('Station');
ylabel('Mole Fraction');
xticks(P7.stations);
plotfixer();