function [Cp, k] = aveSpecHeatAir(T_o, T)
% Calculates the variable specific heat, Cp, for air at a certain temperature in
% Kelvin in J/kgK, Cv in J/kgK and k.
a = 28.11; b = 0.1967e-2; c = 0.4802e-5; d = -1.966e-9;
R = 287; % J/kgK 
mol2kg = 28.9647e-3; %kg/mol

I_1 = ( a.*(T_o - T) + (b./2).*(T_o.^2 - T.^2) + ...
    (c./3).*(T_o.^3 - T.^3) + (d./4).*(T_o.^4 - T.^4));

Cp = I_1 ./ (T_o - T);
Cp = Cp ./ mol2kg;

Cv = Cp - R;

k = Cp / Cv;
end