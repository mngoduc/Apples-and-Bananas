function T_2 = workOut2Temp(m_dot, W_out, T_1)
    LHS = W_out/m_dot;
    
    % initializing stuff
    T_2 = T_1; 
    RHS = dh(T_1, T_2);
    
    dT = 0.001;
    while LHS > RHS
        T_2 = T_2 - dT;
        RHS = dh(T_1, T_2);
    end 
end 