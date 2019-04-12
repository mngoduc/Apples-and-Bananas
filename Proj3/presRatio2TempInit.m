function T_1 = presRatio2TempInit(P_1, P_o1, T_o1)
    
    LHS = P_o1./P_1;
    
    T_1 =  T_o1;
    
    mol2kg = 28.9647e-3;

    a = 28.11/mol2kg; 
    b = 0.1967e-2/mol2kg;
    c = 0.4802e-5/mol2kg; 
    d = -1.966e-9/mol2kg;
    R = 287; %J/kgK 
    
    I = a*log(T_o1./T_1) + b*(T_o1-T_1) + c/2*(T_o1.^2-T_1.^2) ...
        + d/3*(T_o1.^3-T_1.^3);
    RHS = exp(I/R);
    
    dT = 0.01;
    
    while LHS > RHS
        T_1 = T_1 - dT;
        I = a*log(T_o1./T_1) + b*(T_o1-T_1) + c/2*(T_o1.^2-T_1.^2) ...
                + d/3*(T_o1.^3-T_1.^3);
        RHS = exp(I/R);
    end
end