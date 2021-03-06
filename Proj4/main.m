close all; clear all; clc;
%% Initial Set up 
atm2Pa = 101325; 
C2K = 273; 
T_ref = 298; 
R_u = 8.314;                % [J/mol-K]

% given lambda 
G1.L = 2; 
G1.Pm = 1*atm2Pa; 

% Mass of 
M_H2 = 2.016;                               % [g/mol]
M_H2O = 18.01528e-3;                           % [kg/mol]

% LHV, HHV of H2 from Table A-27
C.LHV_bar = 120e3 * M_H2;                                  % [J/g] --> [J/mol]
C.HHV_bar = 141.8e3 * M_H2 ;                               % [J/g] --> [J/mol]

C.coef_H2Ov = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9];
C.coef_H2 = [29.11, -0.1916e-2, 0.4003e-5, -0.8704e-9];
C.coef_O2 = [25.48, 1.520e-2, -0.7155e-5, 1.312e-9];
C.coef_N2 = [28.90, -0.1571e-2, 0.8081e-5, -2.873e-9];

% heat of formation @ STP, N2 
C.h_fo_N2 = 0;                                              % [J/mol]
% entropy of formation @ STP, N2
C.s_o_N2 = 191.61;    
% heat of formation @ STP, O2
C.h_fo_O2 = 0;                                              % [J/mol]
% entropy of formation @ STP, O2
C.s_o_O2 = 205.04;                                          % [J/mol K] 
% heat of formation @ STP, H2
C.h_fo_H2 = 0;                                              % [J/mol]
% entropy of formation @ STP, H2
C.s_o_H2 = 130.68;                                          % [J/mol K]
% H2O vapor 
C.h_fo_H2Ov = -241820;                                      % [J/mol]
C.s_o_H2Ov = 188.83;                                        % [J/mol K]
% H2O liquid
C.h_fo_H2Ol = -285830;                                      % [J/mol]
C.s_o_H2Ol = 69.92;
C.cp_bar_H2Ol = 4.18*M_H2O;

%% Problem 1 
P1.T = [25:1:1000] + C2K;
P1.Psat = T2P_sat(P1.T);

P1.y_max = P1.Psat./G1.Pm; 
P1.N_H2O = 1;       % moles of water = beta + gamma 
P1.y_test = P1.N_H2O/(P1.N_H2O + 0.5*(G1.L-1) + 0.5*G1.L*3.76); %P1.N_H2O/(P1.N_H2O + 0.5*G1.L);

P1.N_a =  0.5*(G1.L-1) + 0.5*G1.L*3.76;
P1.N_react = 1 + 0.5*G1.L*(1+3.76); 
 
[P1.beta, P1.gamma] = vaporLiquidBalance (P1.T, P1.y_test, P1.y_max, ...
                                P1.N_a, P1.N_H2O, 0);

P1.N_prod = P1.beta + 0.5*(G1.L - 1) + 0.5*(G1.L)*(3.76);

% enthalpy of N2 at T
P1.dh_bar_N2 = integral_h (C.coef_N2, P1.T);      % [J/mol]
                         % [J/mol K] 
% N2 on the REACTANTS side
P1.y_N2_react = (.5*G1.L*3.76)./(P1.N_react-1);          % mole fraction of N2 on reactants side
P1.ds_o_N2_react = delta_s (C.coef_N2, P1.T, P1.y_N2_react);  % integral term and log term of entropy 

P1.g_bar_N2react = C.h_fo_N2 + P1.dh_bar_N2 ...
            - P1.T.*(C.s_o_N2 + P1.ds_o_N2_react);

% N2 on the PROD side
P1.y_N2_prod = (.5*G1.L*3.76)./P1.N_prod;          % mole fraction of N2 on PRODUCT side
P1.ds_o_N2_prod = delta_s (C.coef_N2, P1.T, P1.y_N2_prod);  % integral term and log term of entropy 

P1.g_bar_N2prod = C.h_fo_N2 + P1.dh_bar_N2 ...
            - P1.T.*(C.s_o_N2 + P1.ds_o_N2_prod);
        
% enthalpy of O2 at T
P1.dh_bar_O2 = integral_h (C.coef_O2, P1.T);      % [J/mol]

P1.y_O2_react = (.5*G1.L)./(P1.N_react-1);                              % mole fraction of O2 on reactants side
P1.ds_o_O2_react = delta_s (C.coef_O2, P1.T, P1.y_O2_react);  % integral term and log term of entropy 

P1.g_bar_O2react = C.h_fo_O2 + P1.dh_bar_O2 ...
            - P1.T.*(C.s_o_O2 + P1.ds_o_O2_react);

%%%%%% PRODUCT SIDE
P1.y_O2_prod = (0.5*(G1.L-1))./P1.N_prod;                               % mole fraction of O2 on products side
P1.ds_o_O2_prod = delta_s (C.coef_O2, P1.T, P1.y_O2_prod);    % integral term and log term of entropy

P1.g_bar_O2prod = C.h_fo_O2 + P1.dh_bar_O2 ...
            - P1.T.*(C.s_o_O2 + P1.ds_o_O2_prod);

P1.dh_bar_H2 = integral_h (C.coef_H2, P1.T);      % [J/mol]

P1.y_H2 = 1./P1.N_react;
P1.ds_o_H2 = delta_s (C.coef_H2, P1.T, 1);      % integral term and log term of entropy

P1.g_bar_H2 = C.h_fo_H2 + P1.dh_bar_H2 ...
            - P1.T.*(C.s_o_H2 + P1.ds_o_H2);

% H2O vapor 
P1.dh_bar_H2Ov = integral_h (C.coef_H2Ov, P1.T);              % [J/mol]

P1.y_H2Ov = P1.beta./P1.N_prod; 
P1.ds_o_H2Ov = delta_s (C.coef_H2Ov, P1.T, P1.y_H2Ov);        % integral term and log term of entropy
P1.g_bar_H2Ov = C.h_fo_H2Ov + P1.dh_bar_H2Ov ...
                    - P1.T.*(C.s_o_H2Ov + P1.ds_o_H2Ov);

P1.dh_bar_H2Ol = C.cp_bar_H2Ol*(P1.T - T_ref);             % [J/mol]
P1.y_H2Ol = P1.gamma./P1.N_prod;
P1.ds_o_H2Ol = C.cp_bar_H2Ol.*log(P1.T/T_ref); 
P1.g_bar_H2Ol = C.h_fo_H2Ol + P1.dh_bar_H2Ol ...
                         - P1.T.*(C.s_o_H2Ol + P1.ds_o_H2Ol);

% P1.DG = Delta G [J]
P1.DG = P1.beta .* P1.g_bar_H2Ov ...
        + P1.gamma .* P1.g_bar_H2Ol ...
        + 0.5 *(G1.L - 1).* P1.g_bar_O2prod ...
        + (.5*G1.L*3.76).* P1.g_bar_N2prod...
        - P1.g_bar_H2...                                    % N_H2 = 1  
        - 0.5*G1.L*P1.g_bar_O2react...                      % N_O2_react = 1 = 1/2 * 2 
        - (.5*G1.L*3.76).* P1.g_bar_N2react; 

P1.Dh_bar = P1.beta.*(P1.dh_bar_H2Ov + C.h_fo_H2Ov) ...
        + P1.gamma.*(P1.dh_bar_H2Ol + C.h_fo_H2Ol) ...
        + .5*(G1.L-1)*(P1.dh_bar_O2 + C.h_fo_O2) ...
        - 1*(P1.dh_bar_H2 + C.h_fo_H2) ...
        - .5*G1.L*(P1.dh_bar_O2 + C.h_fo_O2);                % [J/mol]

P1.N_fuel = 1;                                              % [mol]   
P1.DH = P1.Dh_bar * P1.N_fuel;                              % [J]
P1.LHV = C.LHV_bar * P1.N_fuel;                            % [J]
P1.HHV = C.HHV_bar * P1.N_fuel;                            % [J]

P1.E_HHV = (-P1.DG)./(P1.HHV);
P1.E_LHV = (-P1.DG)./(P1.LHV);
P1.E_DH = (P1.DG)./(P1.DH);

% Carnot Efficiency  = 1 - TL/TH
P1.E_c = 1 - (25 + C2K)./(P1.T); 

% Plots for part 1
figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
plot(P1.T,P1.E_HHV,'LineWidth',2,'MarkerSize',10);
plot(P1.T,P1.E_LHV,'LineWidth',2,'MarkerSize',10);
plot(P1.T,P1.E_DH,'LineWidth',2,'MarkerSize',10);
plot(P1.T,P1.E_c,'LineWidth',2,'MarkerSize',10);
xlabel('Temperature, [K]','FontSize',18,'FontWeight','bold');
ylabel('First Law Efficiency, \eta_I','FontSize',18,'FontWeight','bold');
P1.leg = legend('HHV','LHV','H_{rxn}', 'Carnot Efficiency', 'Location', 'best');
P1.leg.FontSize = 18; 
set(gca,'FontSize',18);

%% Problem 2

P2.L = linspace(1,10,100);
% mole numbers are unique to each Lambda
P2.N_H2O = 1;       % moles of water = beta + gamma 
P2.y_test = P2.N_H2O./(P2.N_H2O + 0.5*(P2.L-1) + 0.5*P2.L*3.76); %P1.N_H2O/(P1.N_H2O + 0.5*G1.L);

P2.N_a =  0.5*(G1.L-1) + 0.5*G1.L*3.76;
P2.N_react = 1 + 0.5*P2.L*(1+3.76); 

P2.T = [80, 220, 650, 800] + C2K; % 80 + C2K;%

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;
for j = 1:length(P2.T)
    P2.Psat = T2P_sat(P2.T(j));                 % unique to each temp 
    P2.y_max = P2.Psat./G1.Pm;                  % unique to each temp 
[P2.beta, P2.gamma] = vaporLiquidBalance (P2.T(j), P2.y_test, P2.y_max,...
                        P2.N_a, P2.N_H2O, 0);
P2.N_prod = P2.beta + 0.5*(P2.L - 1) + 0.5*(P2.L)*(3.76);

% enthalpy of N2 at T
P2.dh_bar_N2 = integral_h (C.coef_N2, P2.T(j));      % [J/mol]

% N2 on the REACTANTS side
P2.y_N2_react = (.5*G1.L*3.76)./(P2.N_react-1);          % mole fraction of N2 on reactants side
P2.ds_o_N2_react = delta_s (C.coef_N2, P2.T(j), P2.y_N2_react);  % integral term and log term of entropy 

P2.g_bar_N2react(j,:) = C.h_fo_N2 + P2.dh_bar_N2 ...
                            - P2.T(j).*(C.s_o_N2 + P2.ds_o_N2_react);

% N2 on the PROD side
P2.y_N2_prod = (.5*G1.L*3.76)./P2.N_prod;          % mole fraction of N2 on PRODUCT side
P2.ds_o_N2_prod = delta_s (C.coef_N2, P2.T(j), P2.y_N2_prod);  % integral term and log term of entropy 

P2.g_bar_N2prod(j,:) = C.h_fo_N2 + P2.dh_bar_N2 ...
            - P2.T(j).*(C.s_o_N2 + P2.ds_o_N2_prod);

% enthalpy of O2 at T
P2.dh_bar_O2 = integral_h (C.coef_O2, P2.T(j));   % [J/mol]

P2.y_O2_react = (.5*P2.L)./(P2.N_react-1);          % mole fraction of O2 on reactants side
P2.ds_o_O2_react = delta_s (C.coef_O2, P2.T(j), P2.y_O2_react);  % integral term and log term of entropy 

P2.g_bar_O2react(j,:) = C.h_fo_O2 + P2.dh_bar_O2 ...
            - P2.T(j).*(C.s_o_O2 + P2.ds_o_O2_react);
P2.y_O2_prod = (0.5*(P2.L-1))./P2.N_prod;                   % mole fraction of O2 on products side

for k = 1:length(P2.y_O2_prod)
    if P2.y_O2_prod(k) == 0
        P2.g_bar_O2prod(j,:) = 0;
    else
        P2.ds_o_O2_prod = delta_s (C.coef_O2, P2.T(j), P2.y_O2_prod(k));    % integral term and log term of entropy
        P2.g_bar_O2prod(j,k) = C.h_fo_O2 + P2.dh_bar_O2 ...
                        - P2.T(j).*(C.s_o_O2 + P2.ds_o_O2_prod);
    end 
end 

P2.dh_bar_H2 = integral_h (C.coef_H2, P2.T(j));       % [J/mol]

P2.y_H2 = 1./P2.N_react;
P2.ds_o_H2 = delta_s (C.coef_H2, P2.T(j), 1);      % integral term and log term of entropy

P2.g_bar_H2(j,:) = C.h_fo_H2 + P2.dh_bar_H2 ...
                - P2.T(j).*(C.s_o_H2 + P2.ds_o_H2);
% H20 Vapor
P2.dh_bar_H2Ov = integral_h (C.coef_H2Ov, P2.T(j));  % [J/mol]

P2.y_H2Ov = P2.beta./P2.N_prod; 
P2.ds_o_H2Ov = delta_s (C.coef_H2Ov, P2.T(j), P2.y_H2Ov);    % integral term and log term of entropy
P2.g_bar_H2Ov(j,:) = C.h_fo_H2Ov + P2.dh_bar_H2Ov ...
                    - P2.T(j).*(C.s_o_H2Ov + P2.ds_o_H2Ov);

% H2O liquid
P2.dh_bar_H2Ol = C.cp_bar_H2Ol*(P2.T(j) - T_ref);             % [J/mol]
P2.ds_o_H2Ol = C.cp_bar_H2Ol.*log(P2.T(j)/T_ref); 
P2.g_bar_H2Ol(j,:) = ones(size(P2.L)).*(C.h_fo_H2Ol + P2.dh_bar_H2Ol ...
                         - P2.T(j).*(C.s_o_H2Ol + P2.ds_o_H2Ol));

P2.DG(j,:) = P2.beta .* P2.g_bar_H2Ov(j,:) ...
        + P2.gamma .* P2.g_bar_H2Ol(j,:) ...
        + 0.5 *(P2.L - 1).* P2.g_bar_O2prod(j,:) ...
        + (.5*P2.L*3.76).* P2.g_bar_N2prod(j,:)...
        - 1.*P2.g_bar_H2(j,:)...                                   
        - 0.5*P2.L.*P2.g_bar_O2react(j,:)...
        - (.5*P2.L*3.76).* P2.g_bar_N2react(j,:); 

% These are the values that the book gives; Table A-27
P2.LHV_bar = 120e3 * M_H2;
P2.N_fuel = 1;                                              % [mol]   
P2.LHV = P2.LHV_bar * P2.N_fuel; 
P2.E_LHV(j,:) = (-P2.DG(j,:))./(P2.LHV);

plot(P2.L,P2.E_LHV(j,:),'LineWidth',2,'MarkerSize',10);
end
     
xlabel('Excess Air Coefficient, \lambda','FontSize',18,'FontWeight','bold');
ylabel('First Law Efficiency, \eta_I','FontSize',18,'FontWeight','bold');
P2.leg = legend('80^oC', '220^oC', '650^oC', '800^oC', 'Location', 'best');
P2.leg.FontSize = 18; 
set(gca,'FontSize',18);

%% Deliverable 3, Problem 2 - Part 2

P3.P_ref = 1 * atm2Pa; 
P3.Pm = linspace(1,40,100) * atm2Pa;

P3.N_H2O = 1;       % moles of water = beta + gamma 
P3.y_test = P3.N_H2O/(P3.N_H2O + 0.5*(G1.L-1) + 0.5*G1.L*3.76); %P1.N_H2O/(P1.N_H2O + 0.5*G1.L);

P3.N_a =  0.5*(G1.L-1) + 0.5*G1.L*3.76;
P3.N_react = 1 + 0.5*G1.L*(1+3.76); 

P3.T = [80, 220, 650, 800] + C2K;
P3.Psat = T2P_sat(P3.T);                 % unique to each temp 

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on;

for j = 1:length(P3.T)
P3.y_max = P3.Psat(j)./P3.Pm;                  % unique to each temp 

[P3.beta, P3.gamma] = vaporLiquidBalance (P3.T(j), P3.y_test, P3.y_max,...
                            P3.N_a, P3.N_H2O, 0);

P3.N_prod = P3.beta + 0.5*(G1.L - 1) + 0.5*(G1.L)*(3.76);

% enthalpy of N2 at T
P3.dh_bar_N2 = integral_h (C.coef_N2, P3.T(j));      % [J/mol]

% N2 on the REACTANTS side
P3.y_N2_react = (.5*G1.L*3.76)./(P3.N_react-1);          % mole fraction of N2 on reactants side
P3.p_rat_N2_react = P3.y_N2_react.*P3.Pm ./ P3.P_ref;
P3.ds_o_N2_react = delta_s (C.coef_N2, P3.T(j), P3.p_rat_N2_react);  % integral term and log term of entropy 

P3.g_bar_N2react = C.h_fo_N2 + P3.dh_bar_N2 ...
            - P3.T(j).*(C.s_o_N2 + P3.ds_o_N2_react);

% N2 on the PROD side
P3.y_N2_prod = (.5*G1.L*3.76)./P3.N_prod;          % mole fraction of N2 on PRODUCT side
P3.p_rat_N2_prod = P3.y_N2_prod.*P3.Pm ./ P3.P_ref;
P3.ds_o_N2_prod = delta_s (C.coef_N2, P3.T(j), P3.p_rat_N2_prod);  % integral term and log term of entropy 

P3.g_bar_N2prod = C.h_fo_N2 + P3.dh_bar_N2 ...
            - P3.T(j).*(C.s_o_N2 + P3.ds_o_N2_prod);
        
% enthalpy of O2 at T
P3.dh_bar_O2 = integral_h (C.coef_O2, P3.T(j));   % [J/mol]

P3.y_O2_react = (.5*G1.L)./(P3.N_react-1);          % mole fraction of O2 on reactants side
P3.p_rat_O2_react = P3.y_O2_react .* P3.Pm ./ P3.P_ref;
P3.ds_o_O2_react = delta_s (C.coef_O2, P3.T(j), P3.p_rat_O2_react);  % integral term and log term of entropy 

P3.g_bar_O2react(j,:) = C.h_fo_O2 + P3.dh_bar_O2 ...
            - P3.T(j).*(C.s_o_O2 + P3.ds_o_O2_react);

P3.y_O2_prod = (0.5*(G1.L-1))./P3.N_prod;                   % mole fraction of O2 on products side
P3.p_rat_O2_prod = P3.y_O2_prod .* P3.Pm ./ P3.P_ref;

for k = 1:length(P3.y_O2_prod)
    if P3.y_O2_prod(k) == 0
        P3.g_bar_O2prod(j,:) = 0;
    else
        P3.ds_o_O2_prod = delta_s (C.coef_O2, P3.T(j), P3.p_rat_O2_prod(k));    % integral term and log term of entropy
        P3.g_bar_O2prod(j,k) = C.h_fo_O2 + P3.dh_bar_O2 ...
                    - P3.T(j).*(C.s_o_O2 + P3.ds_o_O2_prod);
    end 
end 

P3.dh_bar_H2 = integral_h (C.coef_H2, P3.T(j));       % [J/mol]

P3.y_H2 = 1./P3.N_react;
P3.p_rat_H2 = P3.Pm ./ P3.P_ref;
P3.ds_o_H2 = delta_s (C.coef_H2, P3.T(j), P3.p_rat_H2);      % integral term and log term of entropy

P3.g_bar_H2(j,:) = C.h_fo_H2 + P3.dh_bar_H2 ...
                    - P3.T(j).*(C.s_o_H2 + P3.ds_o_H2);

P3.dh_bar_H2Ov = integral_h (C.coef_H2Ov, P3.T(j));  % [J/mol]

P3.y_H2Ov = P3.beta./P3.N_prod;
P3.p_rat_H2Ov = P3.y_H2Ov .* P3.Pm ./ P3.P_ref;
P3.ds_o_H2Ov = delta_s (C.coef_H2Ov, P3.T(j), P3.p_rat_H2Ov);                 % integral term and log term of entropy
P3.g_bar_H2Ov(j,:) = C.h_fo_H2Ov + P3.dh_bar_H2Ov ...
                    - P3.T(j).*(C.s_o_H2Ov + P3.ds_o_H2Ov);

% H2O liquid
P3.dh_bar_H2Ol = C.cp_bar_H2Ol*(P3.T(j) - T_ref);             % [J/mol]
P3.y_H2Ol = P3.gamma./P3.N_prod;
P3.ds_o_H2Ol = C.cp_bar_H2Ol.*log(P3.T(j)/T_ref); 
P3.g_bar_H2Ol(j,:) = ones(size(G1.L)).*(C.h_fo_H2Ol + P3.dh_bar_H2Ol ...
                         - P3.T(j).*(C.s_o_H2Ol + P3.ds_o_H2Ol));

P3.DG(j,:) = P3.beta .* P3.g_bar_H2Ov(j,:) ...
        + P3.gamma .* P3.g_bar_H2Ol(j,:) ...
        + 0.5 *(G1.L - 1).* P3.g_bar_O2prod(j,:) ...
        + (.5*G1.L*3.76).* P3.g_bar_N2prod...
        - (.5*G1.L*3.76).* P3.g_bar_N2react...
        - 1.*P3.g_bar_H2(j,:)...                                    % N_H2 = 1  
        - 0.5*G1.L.*P3.g_bar_O2react(j,:);  
                     
P3.E_LHV(j,:) = (-P3.DG(j,:))./(P1.LHV);

plot(P3.Pm/atm2Pa,P3.E_LHV(j,:),'LineWidth',2,'MarkerSize',10);
end
     
xlabel('Pressure, [atm]','FontSize',18,'FontWeight','bold');
ylabel('First Law Efficiency, \eta_I','FontSize',18,'FontWeight','bold');
P3.leg = legend('80^oC', '220^oC', '650^oC', '800^oC', 'Location', 'best');
P3.leg.FontSize = 18; 
set(gca,'FontSize',18);
 
%% Deliverable 4, Problem 3 - inlet RH vs temperature required for a saturated outlet air stream

% RH = mv/mg = (mass of vapor)/(mass at sat.)
% RH = Pv/Pg
% beta + gamma = 1 + alpha
P4.Pm = G1.Pm;
P4.L = 2;
P4.T = [25:1:100] + C2K;

P4.Psat = T2P_sat(P4.T);  
P4.y_max = P4.Psat./P4.Pm;

P4.beta = P4.y_max./(1-P4.y_max)*(0.5*(P4.L-1) + 0.5*3.76*P4.L);

for i = 1:length(P4.beta)
    if P4.beta(i) < 1
        P4.alpha(i) = 0;
        P4.gamma(i) = 1-P4.beta(i);
    else
        P4.alpha(i) = P4.beta(i) - 1;
        P4.gamma(i) = 0;
    end
end

P4.Pv = P4.alpha./(0.5*P4.L+0.5*P4.L*3.76+P4.alpha).*P4.Pm;
P4.Pg = P4.Psat;
P4.RH = P4.Pv./P4.Pg;

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
plot(P4.T,P4.RH,'LineWidth',2,'MarkerSize',10);
xlabel('Temperature, [K]','FontSize',18,'FontWeight','bold');
ylabel('Relative Humidity','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18);

%% Deliverable 5

P5.L = 2; 
P5.T = [25:1:100] + C2K;
P5.Pm = 1*atm2Pa; 

figure('units','normalized','outerposition',[0 0 .75 .75]); % for larger plot
hold on
% Line (1): dry hydrogen and dry air,
P5.DG = P1.DG(P1.T <= 100 + C2K);
P5.E_LHV = (-P5.DG)./(P1.LHV);

plot(P5.T,P5.E_LHV,'LineWidth',2,'MarkerSize',10);

%%%%%%%%%% shared constants in parts 2 and 3: %%%%%%%%%%%%%%%
% enthalpy of N2 at T
P5.dh_bar_N2 = integral_h (C.coef_N2, P5.T);        % [J/mol]
% enthalpy of O2 at T
P5.dh_bar_O2 = integral_h (C.coef_O2, P5.T);        % [J/mol]

%%%%%%% H2O %%%%%%%%%
% H2O vapor 
P5.dh_bar_H2Ov = integral_h (C.coef_H2Ov, P5.T);    % [J/mol]
P5.dh_bar_H2Ol = C.cp_bar_H2Ol*(P5.T - T_ref);      % [J/mol]

%%%%%%%% N2 %%%%%%%%%%
P5.dh_bar_H2 = integral_h (C.coef_H2, P5.T);        % [J/mol]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line (2): dry hydrogen and air with 100% relative humidity at cell inlet 
% conditions, 

P5.Psat = T2P_sat(P5.T);
P5.y_max = P5.Psat./P5.Pm;

P5.alpha2 = P5.y_max./(1-P5.y_max).*(0.5*P5.L*(1+ 3.76));

P5.N_a_prod2 =  0.5*(P5.L-1) + 0.5*P5.L*3.76;

P5.N_H2O_prod2 = 1 + P5.alpha2;
P5.y_test_prod = P5.N_H2O_prod2 ./ (P5.N_H2O_prod2 + 0.5*(P5.L-1) + 0.5*P5.L*3.76); 

[P5.beta2, P5.gamma2] = vaporLiquidBalance (P5.T, P5.y_test_prod, ...
                        P5.y_max, P5.N_a_prod2, P5.N_H2O_prod2, P5.alpha2);

P5.N_prod2 = P5.beta2 + 0.5*(P5.L-1) + 0.5*(P5.L)*3.76;
P5.N_react2 = 0.5*(P5.L)*(1+3.76) + P5.alpha2;

% N2 on the REACTANTS side
P5.y_N2_react2 = (.5*G1.L*3.76)./(P5.N_react2);                 % mole fraction of N2 on reactants side
P5.ds_o_N2_react2 = delta_s (C.coef_N2, P5.T, P5.y_N2_react2);    % integral term and log term of entropy 

P5.g_bar_N2react2 = C.h_fo_N2 + P5.dh_bar_N2 ...
            - P5.T.*(C.s_o_N2 + P5.ds_o_N2_react2);

% N2 on the PROD side
P5.y_N2_prod2 = (.5*G1.L*3.76)./P5.N_prod2;                      % mole fraction of N2 on PRODUCT side
P5.ds_o_N2_prod2 = delta_s (C.coef_N2, P5.T, P5.y_N2_prod2);        % integral term and log term of entropy 

P5.g_bar_N2prod2 = C.h_fo_N2 + P5.dh_bar_N2 ...
            - P5.T.*(C.s_o_N2 + P5.ds_o_N2_prod2);

P5.y_O2_react2 = (.5*P5.L)./(P5.N_react2);                              % mole fraction of O2 on reactants side
P5.ds_o_O2_react2 = delta_s (C.coef_O2, P5.T, P5.y_O2_react2);  % integral term and log term of entropy 

P5.g_bar_O2react2 = C.h_fo_O2 + P5.dh_bar_O2 ...
            - P5.T.*(C.s_o_O2 + P5.ds_o_O2_react2);
        
% O2, PROD side
P5.y_O2_prod2 = (0.5*(P5.L-1))./P5.N_prod2;                              % mole fraction of O2 on reactants side
P5.ds_o_O2_prod2 = delta_s (C.coef_O2, P5.T, P5.y_O2_prod2);  % integral term and log term of entropy 

P5.g_bar_O2_prod2 = C.h_fo_O2 + P5.dh_bar_O2 ...
            - P5.T.*(C.s_o_O2 + P5.ds_o_O2_prod2);

% H2, reactants
P5.y_H2 = 1./P5.N_react2;
P5.ds_o_H2 = delta_s (C.coef_H2, P5.T, 1);      % integral term and log term of entropy

P5.g_bar_H2 = C.h_fo_H2 + P5.dh_bar_H2 ...
            - P5.T.*(C.s_o_H2 + P5.ds_o_H2);

% H2O vapor reactants
P5.y_H2Ov_react2 = P5.alpha2./P5.N_react2;
P5.ds_o_H2Ov_react2 = delta_s (C.coef_H2Ov, P5.T, P5.y_H2Ov_react2);        % integral term and log term of entropy
P5.g_bar_H2Ov_react2 = C.h_fo_H2Ov + P5.dh_bar_H2Ov ...
                    - P5.T.*(C.s_o_H2Ov + P5.ds_o_H2Ov_react2);
        
% H2O vapor products    
P5.y_H2Ov_prod2 = P5.beta2./P5.N_prod2; 
P5.ds_o_H2Ov_prod2 = delta_s (C.coef_H2Ov, P5.T, P5.y_H2Ov_prod2);        % integral term and log term of entropy
P5.g_bar_H2Ov_prod2 = C.h_fo_H2Ov + P5.dh_bar_H2Ov ...
                    - P5.T.*(C.s_o_H2Ov + P5.ds_o_H2Ov_prod2);

% H2O Liquid products
P5.ds_o_H2Ol_prod2 = C.cp_bar_H2Ol.*log(P5.T/T_ref); 
P5.g_bar_H2Ol_prod2 = C.h_fo_H2Ol + P5.dh_bar_H2Ol ...
                         - P5.T.*(C.s_o_H2Ol + P5.ds_o_H2Ol_prod2);

% P5.DG = Delta G [J]
P5.DG2 = P5.beta2 .* P5.g_bar_H2Ov_prod2 ...
        + P5.gamma2 .* P5.g_bar_H2Ol_prod2 ...
        + 0.5 *(P5.L - 1).* P5.g_bar_O2_prod2 ...
        + (.5*P5.L*3.76).* P5.g_bar_N2prod2...
        - P5.g_bar_H2...                                    % N_H2 = 1  
        - 0.5*P5.L*P5.g_bar_O2react2...                      % N_O2_react = 1 = 1/2 * 2 
        - (.5*P5.L*3.76).* P5.g_bar_N2react2...
        - P5.alpha2.*P5.g_bar_H2Ov_react2; 
    
P5.P_H2O_v = P5.alpha2./(P5.alpha2+.5.*P5.L.*(1+3.76)).*P5.Pm;
P5.RH = P5.P_H2O_v./P5.Psat;
                     
P5.E_LHV2 = (-P5.DG2)./(P1.LHV);
plot(P5.T,P5.E_LHV2, ':','LineWidth',2,'MarkerSize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line (3): dry hydrogen and air with just enough input humidity added to 
% maintain saturation conditions at the exit.
P5.alpha3 = P4.alpha; 
P5.beta3 = P4.beta;
P5.gamma = P4.gamma; 

P5.N_react3 = 0.5*P5.L*(1+3.76) + P5.alpha3; 
P5.N_prod3 = P5.beta3 + 0.5*(P5.L - 1) + 0.5*(P5.L)*(3.76);

% entropy of N2 at T, P
% N2 on the REACTANTS side
P5.y_N2_react = (.5*P5.L*3.76)./(P5.N_react3);          % mole fraction of N2 on reactants side
P5.ds_o_N2_react = delta_s (C.coef_N2, P5.T, P5.y_N2_react);  % integral term and log term of entropy 

P5.g_bar_N2react = C.h_fo_N2 + P5.dh_bar_N2 ...
            - P5.T.*(C.s_o_N2 + P5.ds_o_N2_react);

% N2 on the PROD side
P5.y_N2_prod = (.5*G1.L*3.76)./P5.N_prod3;          % mole fraction of N2 on PRODUCT side
P5.ds_o_N2_prod = delta_s (C.coef_N2, P5.T, P5.y_N2_prod);  % integral term and log term of entropy 

P5.g_bar_N2prod = C.h_fo_N2 + P5.dh_bar_N2 ...
            - P5.T.*(C.s_o_N2 + P5.ds_o_N2_prod);

%%%%%% REACTANT SIDE, O2 %%%%%%
P5.y_O2_react = (.5*G1.L)./(P5.N_react3);                              % mole fraction of O2 on reactants side
P5.ds_o_O2_react = delta_s (C.coef_O2, P5.T, P5.y_O2_react);  % integral term and log term of entropy 

P5.g_bar_O2react = C.h_fo_O2 + P5.dh_bar_O2 ...
            - P5.T.*(C.s_o_O2 + P5.ds_o_O2_react);

%%%%%% PRODUCT SIDE, O2 %%%%%%
P5.y_O2_prod = (0.5*(G1.L-1))./P5.N_prod3;                               % mole fraction of O2 on products side
P5.ds_o_O2_prod = delta_s (C.coef_O2, P5.T, P5.y_O2_prod);    % integral term and log term of entropy

P5.g_bar_O2prod = C.h_fo_O2 + P5.dh_bar_O2 ...
            - P5.T.*(C.s_o_O2 + P5.ds_o_O2_prod);

%%%%% REACTANT SIDE, H2 %%%%%%%%
P5.ds_o_H2 = delta_s (C.coef_H2, P5.T, 1);      % integral term and log term of entropy
P5.g_bar_H2 = C.h_fo_H2 + P5.dh_bar_H2 ...
            - P5.T.*(C.s_o_H2 + P5.ds_o_H2);

%%%%%%% PRODUCT SIDE, H2O vapor %%%%%%%%%
P5.y_H2Ov_prod = P5.beta3./P5.N_prod3; 
P5.ds_o_H2Ov_prod = delta_s (C.coef_H2Ov, P5.T, P5.y_H2Ov_prod);        % integral term and log term of entropy
P5.g_bar_H2Ov_prod = C.h_fo_H2Ov + P5.dh_bar_H2Ov ...
                    - P5.T.*(C.s_o_H2Ov + P5.ds_o_H2Ov_prod);

% H2O liquid
P5.y_H2Ol = P5.gamma./P5.N_prod3;
P5.ds_o_H2Ol = C.cp_bar_H2Ol.*log(P5.T/T_ref); 
P5.g_bar_H2Ol = C.h_fo_H2Ol + P5.dh_bar_H2Ol ...
                         - P5.T.*(C.s_o_H2Ol + P5.ds_o_H2Ol);

%%%%%%%%% REACTANT SIDE H2O %%%%%%%%%%
P5.y_H2Ov_r = P5.alpha3./P5.N_react3; 
for i = 1:length(P5.y_H2Ov_r)
    if P5.y_H2Ov_r(i) == 0
        P5.ds_o_H2Ov_r(i) = 0; 
    else
        P5.ds_o_H2Ov_r(i) = delta_s (C.coef_H2Ov, P5.T(i), P5.y_H2Ov_r(i));        % integral term and log term of entropy
    end
end

P5.g_bar_H2Ov_react = C.h_fo_H2Ov + P5.dh_bar_H2Ov ...
                    - P5.T.*(C.s_o_H2Ov  + P5.ds_o_H2Ov_r);

% P5.DG = Delta G [J]
P5.DG_sat_out = P5.beta3 .* P5.g_bar_H2Ov_prod ...
        + P5.gamma .* P5.g_bar_H2Ol ...
        + 0.5 *(G1.L - 1).* P5.g_bar_O2prod ...
        + (.5*G1.L*3.76).* P5.g_bar_N2prod...
        - P5.g_bar_H2...                                    % N_H2 = 1  
        - 0.5*G1.L*P5.g_bar_O2react...                      % N_O2_react = 1 = 1/2 * 2 
        - (.5*G1.L*3.76).* P5.g_bar_N2react...
        - P5.alpha3.* P5.g_bar_H2Ov_react; 

P5.E_sat_out = -(P5.DG_sat_out)./(P1.LHV);
plot(P5.T,P5.E_sat_out, ':','LineWidth',2,'MarkerSize',10);

xlabel('Temperature, [K]','FontSize',18,'FontWeight','bold');
ylabel('First Law Efficiency, \eta_I','FontSize',18,'FontWeight','bold');
P5.leg = legend('Dry H_2, Dry Air','Saturated Inlet Air','Saturated Exit Air','Location', 'best');
P5.leg.FontSize = 18; 
set(gca,'FontSize',18);
