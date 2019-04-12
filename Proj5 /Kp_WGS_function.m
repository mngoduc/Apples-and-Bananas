function Kp_WGS = Kp_WGS_function(T)

% WGS: H2O + CO -> CO2 + H2
R_u = 8.314;                      % [J/ mol K]
coef_H2Ov = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9];      % should be [J/mol]
coef_H2 = [29.11, -0.1916e-2, 0.4003e-5, -0.8704e-9];
coef_CO = [28.16, 0.1675e-2, 0.5372e-5, -2.222e-9];
coef_CO2 = [22.26, 5.981e-2, -3.501e-5, 7.469e-9];

% heat of formation @ STP, H2
h_fo_H2 = 0;                                                % [J/mol]
% entropy of formation @ STP, H2
s_o_H2 = 130.68;                                            % [J/mol K]
% heat of formation @ STP, H2O vapor 
h_fo_H2Ov = -241820;                                        % [J/mol]
% entropy of formation @ STP, H2O vapor
s_o_H2Ov = 188.83;                                          % [J/mol K]
% heat of formation @ STP, CO
h_fo_CO = -110530;                                          % [J/mol]
% entropy of formation @ STP, CO
s_o_CO = 197.65;                                            % [J/mol K]
% heat of formation @ STP, CO2
h_fo_CO2 = -393520;                                         % [J/mol]
% entropy of formation @ STP, CO2
s_o_CO2 = 213.8;                                            % [J/mol K]

% H2Ov WGS REACTANT 
N_H2Ov_WGS = 1;                                             % number of mols
% y_H2Ov_WGS =  N_H2Ov_WGS ./ N_WGS_R;                            
% mole fraction; assume Pm = 1 atm 
dh_bar_H2Ov_WGS = integral_h ( coef_H2Ov,  T);              % [J/mol]
ds_o_H2Ov_WGS = delta_s ( coef_H2Ov,  T, 1);                % y_H2Ov_WGS);          
% integral term and log term of entropy 
g_bar_H2Ov_WGS =  h_fo_H2Ov +  dh_bar_H2Ov_WGS ...
                -  T.*( s_o_H2Ov +  ds_o_H2Ov_WGS);

% CO WGS REACTANT 
N_CO_WGS = 1;
%  y_CO_WGS =  N_CO_WGS ./ N_WGS_R;                            
% mole fraction; assume Pm = 1 atm 
dh_bar_CO_WGS = integral_h ( coef_CO,  T);
ds_o_CO_WGS = delta_s ( coef_CO,  T, 1);% y_CO_WGS);          
% integral term and log term of entropy 
g_bar_CO_WGS =  h_fo_CO +  dh_bar_CO_WGS ...
                -  T.*( s_o_CO+  ds_o_CO_WGS);              % [J/mol]

% CO2 WGS PRODUCT
N_CO2 = 1;
%  y_CO2 =  N_CO2 ./ N_WGS_P;                            
% mole fraction; assume Pm = 1 atm 
dh_bar_CO2 = integral_h ( coef_CO2,  T);                    % [J/mol]
ds_o_CO2 = delta_s ( coef_CO2,  T, 1);% y_CO2);          
% integral term and log term of entropy 
g_bar_CO2 =  h_fo_CO2 +  dh_bar_CO2 ...
                -  T.*( s_o_CO2 +  ds_o_CO2);               % [J/mol]

% H2 WGS PRODUCT
N_H2_WGS = 1;
%  y_H2_WGS =  N_H2_WGS ./ N_WGS_P;                            
% mole fraction; assume Pm = 1 atm 
dh_bar_H2_WGS = integral_h ( coef_H2,  T);                  
ds_o_H2_WGS = delta_s ( coef_H2,  T, 1);% y_H2_WGS);          
% integral term and log term of entropy 
g_bar_H2_WGS =  h_fo_H2 +  dh_bar_H2_WGS ...
                -  T.*( s_o_H2 +  ds_o_H2_WGS);             % [J/mol]

% GATHER TERMS FOR WGS 
DGT_WGS =  N_H2_WGS *  g_bar_H2_WGS ...
            +  N_CO2 *  g_bar_CO2 ...
            -  N_CO_WGS *  g_bar_CO_WGS ...
            -  N_H2Ov_WGS *  g_bar_H2Ov_WGS;                % [J] per reaction     

Kp_WGS = exp(- DGT_WGS ./ ( R_u.* T));

end