function Psat = T2P_sat(T)
% temp in Kelvin 
Psat = exp(-1.2914e8.*T.^-3 + 8.2048e5.*T.^-2 - 6522.8./T + 25.5887);
end 