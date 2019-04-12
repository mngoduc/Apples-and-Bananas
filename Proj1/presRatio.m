function [presRatio] = presRatio(k, mach)
% Calculates the pressure ratio based on the K and mach number of a flow
exponent = k ./ (k-1);
C = (k-1)./2;

presRatio = (1 + C.* mach.^2) .^ exponent;
end