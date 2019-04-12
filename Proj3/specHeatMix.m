function [totCp, k] = specHeatMix(T, molFrac)
% Calculates the specific heat, Cp, for air at a certain temperature in
% Kelvin in J/kgK, Cv in J/kgK and k.
% Cp = a + b.*T + c.*T.^2 + d.*T.^3;

    % format = [a, b, c, d]   
    % cp coefficients water(l); [kJ/kmol K]
    coef_w = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9];
    Cp(1,:) = coef_w(1) + coef_w(2)*T + coef_w(3).*T.^2 + coef_w(4).*T.^3;
    
    % cp coefficients carbon dioxide; [kJ/kmol K]
    coef_CO2 = [22.26, 5.981e-2, -3.501e-5, 7.469e-9];
    Cp(2,:) = coef_CO2(1) + coef_CO2(2)*T + coef_CO2(3).*T.^2 + coef_CO2(4).*T.^3;
    
    % cp coefficients N2; [kJ/kmol K]
    coef_N2 = [28.90, -0.1571e-2, 0.8081e-5, -2.873e-9];
    Cp(3,:) = coef_N2(1) + coef_N2(2)*T + coef_N2(3).*T.^2 + coef_N2(4).*T.^3;
    
    % cp coefficients O2; [kJ/kmol K]
    coef_O2 = [25.48, 1.520e-2, -0.7155e-5, 1.312e-9];
    Cp(4,:) = coef_O2(1) + coef_O2(2)*T + coef_O2(3).*T.^2 + coef_O2(4).*T.^3;
    
    
    % Calculating actual CP using mole fraction data [kJ/kmol K]
    for i = 1:length(T)
        totCp(i) = sum(molFrac(:,i) .* Cp(:,i));
    end
    
    R_u = 8.314; % [kJ/kmol K]
   
    totCv = totCp - R_u; % [kJ/kmol K]

    k = totCp ./ totCv;
    totCp = totCp / 0.17; % [J/mol K] * [1 mol/0.170 kg] = [J / kg K]
end

