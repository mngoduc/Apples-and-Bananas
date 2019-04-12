function [tempRatio] = tempRatio(k, mach)
% Calculates the pressure ratio based on the K and mach number of a flow
C = (k-1)./2;

tempRatio = 1 + C.* mach.^2;
end