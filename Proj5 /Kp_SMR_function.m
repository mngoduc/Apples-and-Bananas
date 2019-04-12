function Kp_SMR = Kp_SMR_function(T)

coef_H2Ov = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9];      % should be [J/mol]
coef_H2 = [29.11, -0.1916e-2, 0.4003e-5, -0.8704e-9];
coef_CO = [28.16, 0.1675e-2, 0.5372e-5, -2.222e-9];
coef_CH4 = [19.89, 5.024e-2, 1.269e-5, -11.01e-9];

% heat of formation @ STP, H2
h_fo_H2 = 0;                                              % [J/mol]
% entropy of formation @ STP, H2
s_o_H2 = 130.68;                                          % [J/mol K]
% H2O vapor 
h_fo_H2Ov = -241820;                                      % [J/mol]
s_o_H2Ov = 188.83;                                        % [J/mol K]
% CO (SMR, WGS)
h_fo_CO = -110530;                                        % [J/mol]
s_o_CO = 197.65;                                          % [J/mol K]
% CH4 (WGS)
h_fo_CH4 = -74850;                                       % [J/mol]
s_o_CH4 = 186.16;                                         % [J/mol K]

% SMR : CH4 + H2O -> CO + 3 H2
R_u = 8.314;

% CH4,  SMR REACTANT
N_CH4 = 1;
%  y_CH4 =  N_CH4 ./ N_SMR_R;                            % mole fraction; assume Pm = 1 atm 
dh_bar_CH4 = integral_h ( coef_CH4,  T);                   % [J/mol]
ds_o_CH4 = delta_s ( coef_CH4,  T, 1);% y_CH4);          % integral term and log term of entropy 
g_bar_CH4 =  h_fo_CH4 +  dh_bar_CH4 ...
                -  T.*( s_o_CH4 +  ds_o_CH4);

% Steam, H2Ov, SMR REACTANT 
N_H2Ov_SMR = 1;
%  y_H2Ov_SMR =  N_H2Ov_SMR ./ N_SMR_R;                         % mole fraction; assume Pm = 1 atm 
dh_bar_H2Ov_SMR = integral_h ( coef_H2Ov,  T);                % [J/mol]
ds_o_H2Ov_SMR = delta_s ( coef_H2Ov,  T, 1);% y_H2Ov_SMR);  % integral term and log term of entropy 
g_bar_H2Ov_SMR =  h_fo_H2Ov +  dh_bar_H2Ov_SMR ...
                -  T.*( s_o_H2Ov +  ds_o_H2Ov_SMR);

% CO, SMR PRODUCT
N_CO_SMR = 1;
%  y_CO_SMR =  N_CO_SMR ./ N_SMR_P;                            % mole fraction; assume Pm = 1 atm 
dh_bar_CO_SMR = integral_h ( coef_CO,  T);                   % [J/mol]
ds_o_CO_SMR = delta_s ( coef_CO,  T, 1);% y_CO_SMR);          % integral term and log term of entropy 
g_bar_CO_SMR =  h_fo_CO +  dh_bar_CO_SMR ...
                -  T.*( s_o_CO+  ds_o_CO_SMR);

% H2, SMR PRODUCT 
N_H2_SMR = 3;
%  y_H2_SMR =  N_H2_SMR ./ N_SMR_P;                            % mole fraction; assume Pm = 1 atm 
dh_bar_H2_SMR = integral_h ( coef_H2,  T);                   % [J/mol]
ds_o_H2_SMR = delta_s ( coef_H2,  T, 1);% y_H2_SMR);          % integral term and log term of entropy 
g_bar_H2_SMR =  h_fo_H2 +  dh_bar_H2_SMR ...
                -  T.*( s_o_H2 +  ds_o_H2_SMR);

% GATHER TERMS FOR SMR 
DGT_SMR =  N_H2_SMR *  g_bar_H2_SMR...
            +  N_CO_SMR *  g_bar_CO_SMR ...
            -  N_H2Ov_SMR *  g_bar_H2Ov_SMR...
            -  N_CH4 *  g_bar_CH4;

Kp_SMR = exp(-DGT_SMR ./ (R_u* T));

end