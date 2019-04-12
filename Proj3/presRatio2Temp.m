function T_o2 = presRatio2Temp(P_o1, P_o2, T_o1)
    LHS = P_o2./P_o1;
    
    T_o2 =  T_o1;
    
    mol2kg = 28.9647e-3;

    a = 28.11/mol2kg; 
    b = 0.1967e-2/mol2kg;
    c = 0.4802e-5/mol2kg; 
    d = -1.966e-9/mol2kg;
    R = 287; %J/kgK 
    
    I = a*log(T_o2./T_o1) + b*(T_o2-T_o1) + c/2*(T_o2.^2-T_o1.^2) ...
        + d/3*(T_o2.^3-T_o1.^3);
    RHS = exp(I/R);
    
    dT = 0.0001;
    
    while abs(LHS - RHS) > 0.01
        T_o2 = T_o2 + dT;
        I = a*log(T_o2./T_o1) + b*(T_o2-T_o1) + c/2*(T_o2.^2-T_o1.^2) ...
                + d/3*(T_o2.^3-T_o1.^3);
        RHS = exp(I/R);
    end
end