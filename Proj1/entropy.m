function s = entropy(T_Mat, P_Mat)
    % Universal Gas Constant 
    R = 287; % [J/kg K]
    % Polynomial for Cp in CengelBoles are in units of kJ/kmolK, we need
    % J/kgK
    mol2kg = 28.9647e-3;
    
    % From CengelBoles Appendix A-2, Coefficients of Cp
    a = 28.11/mol2kg; 
    b = 0.1967e-2/mol2kg;
    c = 0.4802e-5/mol2kg; 
    d = -1.966e-9/mol2kg;
    
    s = [0,0];
    for i = 2:length(T_Mat)
        T1 = T_Mat(i-1,:);
        T2 = T_Mat(i,:);
        P1 = P_Mat(i-1,:);
        P2 = P_Mat(i,:);
        
        % Term 1 is the integral term in the 
        term1 = a * log(T2./T1) + b * (T2-T1) + c/2 * (T2.^2-T1.^2) ...
            + d/3 * (T2.^3-T1.^3);
        term2 = R.*log(P2./P1);
        ds = term1 - term2; 
        s(i,:) = s(i-1,:) + ds; 
    end 
end 