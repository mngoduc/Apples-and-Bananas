function [Cp, k] = specHeatAir(T)
% Calculates the specific heat, Cp, for air at a certain temperature in
% Kelvin in J/kgK, Cv in J/kgK and k.
a = 28.11; 
b = 0.1967e-2; 
c = 0.4802e-5; 
d = -1.966e-9;
R = 287; %J/kgK 
mol2kg = 28.9647e-3;

Cp = a + b.*T + c.*T.^2 + d.*T.^3;
Cp = Cp ./ mol2kg;

Cv = Cp - R;

k = Cp./ Cv;
end

