function Tf = heat2Tout(Ti, Q, mDot)
    LHS = Q/mDot;
    
    % initializing stuff
    Tf = Ti; 
    RHS = dh(Tf, Ti);
    
    dT = 0.001;
    while LHS > RHS
        Tf = Tf + dT;
        RHS = dh(Tf, Ti);
    end 
end 