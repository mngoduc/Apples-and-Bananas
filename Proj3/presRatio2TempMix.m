function T_o2 = presRatio2TempMix(P_o1, P_o2, T_o1, molFrac)
    
    M_O2 = 32e-3;
    M_N2 = 28.013e-3;
	M_C = 12.01e-3;
	M_H = 1.01e-3;
    M_H2O = 2*M_H + 0.5*M_O2;
    M_CO2 = M_C + M_O2;
    
    M = [M_CO2, M_H2O, M_N2, M_O2]';
   
    % moleFrac = row/column of mol frac 
    M_Mix= sum(molFrac .* M);
    
    
    % cp coefficients carbon dioxide; [kJ/kmol K]
    coef_CO2 = [22.26, 5.981e-2, -3.501e-5, 7.469e-9];

    coef_H2O = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9];

    % cp coefficients N2; [kJ/kmol K]
    coef_N2 = [28.90, -0.1571e-2, 0.8081e-5, -2.873e-9];
    
    % cp coefficients O2; [kJ/kmol K]
    coef_O2 = [25.48, 1.520e-2, -0.7155e-5, 1.312e-9];

    coef_mix = molFrac(1)*coef_CO2 ...
                + molFrac(2)*coef_H2O...
                + molFrac(3)*coef_N2...
                + molFrac(4)*coef_O2;
    coef_mix = coef_mix/M_Mix;
    
    a = coef_mix(1); 
    b = coef_mix(2); 
    c = coef_mix(3); 
    d = coef_mix(4);
    
    LHS = P_o2./P_o1;
    
    T_o2 =  T_o1;
    R_u = 8.314; 
    R = R_u/M_Mix; %J/kgK 
    
    I = a*log(T_o2./T_o1) + b*(T_o2-T_o1) + c/2*(T_o2.^2-T_o1.^2) ...
        + d/3*(T_o2.^3-T_o1.^3);
    RHS = exp(I/R);
    
    dT = 0.01;
    
    while LHS < RHS
        T_o2 = T_o2 - dT;
            I = a*log(T_o2./T_o1) + b*(T_o2-T_o1) + c/2*(T_o2.^2-T_o1.^2) ...
        + d/3*(T_o2.^3-T_o1.^3);
        RHS = exp(I/R);
    end
    
end