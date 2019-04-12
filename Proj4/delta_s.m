function ds = delta_s (coef, T2, p_ratio)
a = coef(1); 
b = coef(2); 
c = coef(3);
d = coef(4);

T1 = 298;               % [K]

R_u = 8.314;

% the integral term in the entropy equation 
I = a.*log(T2./T1)...
    + b.*(T2 - T1)...
    + c.*(T2.^2 - T1.^2)...
    + d.*(T2.^3 - T1.^3);

% in problems 1 and 2
% log term in equation, mole fraction = pressure ratio
% since P_i = y_i * P_m = y_i * P_o b/c P_m = P_o
% R ln(P_i/P_o)
L = R_u.*log(p_ratio);

ds = I - L;
end 