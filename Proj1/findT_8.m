function T_8 = findT_8 (V_8, T_8s, T_o8)
    % initialization 
    T_8 = T_8s;
    LHS = 0.5*V_8^2; 
    RHS = dh(T_o8, T_8); 
    
    dT = 0.01;
    
    while LHS < RHS
        T_8 = T_8 + dT;
        RHS = dh(T_o8, T_8);
    end 
end 