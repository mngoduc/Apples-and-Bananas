function [h] = dh(T2, T1)
    % Calculates the variable specific heat, Cp, for air at a certain 
    % temperature in Kelvin in J/kgK, Cv in J/kgK and k.
    a = 28.11; b = 0.1967e-2; c = 0.4802e-5; d = -1.966e-9;
    mol2kg = 28.9647e-3; %kg/mol
    %R = 287; %J/kg-K

    h = a*(T2 - T1) ...
        + b/2*(T2.^2 - T1.^2)...
        + c/3*(T2.^3 - T1.^3)...
        + d/4*(T2.^4 - T1.^4); %[J/mol]
    
    h = h/mol2kg;  %[J/kg], unit conversions
end