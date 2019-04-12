function [h_mix] = dhMix(molFrac, T2, T1)
    % Calculates the variable specific heat, Cp, for air at a certain 
    % temperature in Kelvin in J/kgK, Cv in J/kgK and k.
    % molFrac = [CO2; H2O; N2; O2]
    % format = [a, b, c, d]   
    % cp coefficients water(l); [kJ/kmol K]
    
    
    M_O2 = 32e-3;
    M_N2 = 28.013e-3;
	M_C = 12.01e-3;
	M_H = 1.01e-3;
    M_H2O = 2*M_H + 0.5*M_O2;
    M_CO2 = M_C + M_O2;
    M = [M_CO2, M_H2O, M_N2, M_O2]';
   
    for i = 1:length(T1)
        M_Mix(i) = sum(molFrac(:,i) .* M);
    end
   
    % cp coefficients carbon dioxide; [kJ/kmol K]
    coef_CO2 = [22.26, 5.981e-2, -3.501e-5, 7.469e-9];
    h_bar(1,:) = coef_CO2(1).*(T2 - T1) ...
        + coef_CO2(2)./2*(T2.^2 - T1.^2)...
        + coef_CO2(3)./3*(T2.^3 - T1.^3)...
        + coef_CO2(4)./4*(T2.^4 - T1.^4); %[J/mol]
    %h(1,:) = h_bar_CO2/M_CO2; 
    
    coef_H2O = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9];
    h_bar(2,:) = coef_H2O(1).*(T2 - T1) ...
        + coef_H2O(2)./2*(T2.^2 - T1.^2)...
        + coef_H2O(3)./3*(T2.^3 - T1.^3)...
        + coef_H2O(4)./4*(T2.^4 - T1.^4); %[J/mol]
    % h(2,:) = h_bar_H2O/M_H2O;
    
    % cp coefficients N2; [kJ/kmol K]
    coef_N2 = [28.90, -0.1571e-2, 0.8081e-5, -2.873e-9];
    h_bar(3,:) = coef_N2(1).*(T2 - T1) ...
        + coef_N2(2)./2*(T2.^2 - T1.^2)...
        + coef_N2(3)./3*(T2.^3 - T1.^3)...
        + coef_N2(4)./4*(T2.^4 - T1.^4); %[J/mol]
    % h(3,:) = h_bar_N2/M_N2;
    
    % cp coefficients O2; [kJ/kmol K]
    coef_O2 = [25.48, 1.520e-2, -0.7155e-5, 1.312e-9];
    h_bar(4,:) = coef_O2(1).*(T2 - T1) ...
        + coef_O2(2)./2*(T2.^2 - T1.^2)...
        + coef_O2(3)./3*(T2.^3 - T1.^3)...
        + coef_O2(4)./4*(T2.^4 - T1.^4); %[J/mol]
    %h(4,:) = h_bar_O2/M_O2;

    for i = 1:length(T1)
        h_bar_mix(i) = sum(molFrac(:,i) .* h_bar(:,i));
    end
    h_mix = h_bar_mix./M_Mix;
end