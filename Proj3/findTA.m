function T_A = findTA(H_bar_jetA, H_bar_CO2, ...
                                            H_bar_H2O, Phi, T_guess)
%ONLY ACCEPTS SCALARS
T_ref = 298; %K

coef_CO2 = [22.26, 5.981e-2, -3.501e-5, 7.469e-9];
coef_H2O = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9];
coef_N2 = [28.90, -0.1571e-2, 0.8081e-5, -2.873e-9];
coef_O2 = [25.48, 1.520e-2, -0.7155e-5, 1.312e-9];

c = 11.1;
d = 17.85*3.76;
a = 17.85;
b = 12.3;

T1 = T_ref;

RHS = (Phi * H_bar_jetA);

% stuff for the denominator
dh_CO2 = @(T2) coef_CO2(1).*(T2 - T1) ...
    + coef_CO2(2)./2*(T2.^2 - T1.^2)...
    + coef_CO2(3)./3*(T2.^3 - T1.^3)...
    + coef_CO2(4)./4*(T2.^4 - T1.^4); %[J/mol]

dh_H2O = @(T2) coef_H2O(1).*(T2 - T1) ...
    + coef_H2O(2)./2*(T2.^2 - T1.^2)...
    + coef_H2O(3)./3*(T2.^3 - T1.^3)...
    + coef_H2O(4)./4*(T2.^4 - T1.^4); %[J/mol]    

dh_N2 = @(T2)  coef_N2(1).*(T2 - T1) ...
    + coef_N2(2)./2*(T2.^2 - T1.^2)...
    + coef_N2(3)./3*(T2.^3 - T1.^3)...
    + coef_N2(4)./4*(T2.^4 - T1.^4); %[J/mol]

dh_O2 = @(T2)  coef_O2(1).*(T2 - T1) ...
    + coef_O2(2)./2*(T2.^2 - T1.^2)...
    + coef_O2(3)./3*(T2.^3 - T1.^3)...
    + coef_O2(4)./4*(T2.^4 - T1.^4); %[J/mol]  

LHS = @(T2) (Phi*b*(H_bar_CO2 + dh_CO2(T2)) + c*Phi*(H_bar_H2O + dh_H2O(T2)) ...
        + d * dh_N2(T2) + (a*(1-Phi)) * dh_O2(T2));
    
fun = @(T2) RHS - LHS(T2);

% options = optimset('Display','iter');
% options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
T_A = fzero(fun, T_guess);


end
 