function deltaS = deltaEntropy(T1, T2, P1, P2)
    cp_avg = aveSpecHeatAir(T2, T1);
    R = 287; % [J/kg K]
    deltaS = cp_avg*log(T2./T1) - R*log(P2./P1);
end 