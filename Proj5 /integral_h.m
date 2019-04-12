function delta_h_bar = integral_h (coef, T2)

% coef = [a,b,c,d]
a = coef(1); 
b = coef(2);
c = coef(3);
d = coef(4);

T1 = 298;                               %[K]

delta_h_bar = a.*(T2 - T1) ...
            + b/2.*(T2.^2 - T1^2)...
            + c/3.*(T2.^3 - T1^3)...
            + d/4.*(T2.^4 - T1^4);       %[J/mol]
end 