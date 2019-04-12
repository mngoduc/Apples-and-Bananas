function [T_a] = findT_a(H_bar_jetA, H_bar_CO2, H_bar_H2O, Phi)
%ONLY ACCEPTS SCALARS
T_ref = 298; %K
T_a = 350;%T_ref; % first guess for adiabatic temp
LHS = T_a;

T_react = 300; %K

c = 11.1;
d = 17.85*3.76;
a = 17.85;

T1 = T_ref;
T2 = T_a;

coef_CO2 = [22.26, 5.981e-2, -3.501e-5, 7.469e-9];
    dh_CO2 = coef_CO2(1).*(T2 - T1) ...
        + coef_CO2(2)./2*(T2.^2 - T1.^2)...
        + coef_CO2(3)./3*(T2.^3 - T1.^3)...
        + coef_CO2(4)./4*(T2.^4 - T1.^4); %[J/mol]
    
coef_w = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9];
    dh_w = coef_w(1).*(T2 - T1) ...
        + coef_w(2)./2*(T2.^2 - T1.^2)...
        + coef_w(3)./3*(T2.^3 - T1.^3)...
        + coef_w(4)./4*(T2.^4 - T1.^4); %[J/mol]
    
     coef_N2 = [28.90, -0.1571e-2, 0.8081e-5, -2.873e-9];
    dh_N2 = coef_N2(1).*(T2 - T1) ...
        + coef_N2(2)./2*(T2.^2 - T1.^2)...
        + coef_N2(3)./3*(T2.^3 - T1.^3)...
        + coef_N2(4)./4*(T2.^4 - T1.^4); %[J/mol]
    
    coef_O2 = [25.48, 1.520e-2, -0.7155e-5, 1.312e-9];
     dh_O2 = coef_O2(1).*(T2 - T1) ...
        + coef_O2(2)./2*(T2.^2 - T1.^2)...
        + coef_O2(3)./3*(T2.^3 - T1.^3)...
        + coef_O2(4)./4*(T2.^4 - T1.^4); %[J/mol]
    
    denom = Phi*(H_bar_CO2 + dh_CO2) + c*Phi*(H_bar_H2O + dh_w) ...
        + d* dh_N2 + (a - Phi*(2+c))*dh_O2;
    num = Phi* H_bar_jetA*T_react;
    
    RHS = num/denom;
    
    dT = 0.00001;
    
    while abs(LHS - RHS) > 0.001
        T_a = T_a +dT;
        
        T2 = T_a;
        
        LHS = T_a;
        
    dh_CO2 = coef_CO2(1).*(T2 - T1) ...
        + coef_CO2(2)./2*(T2.^2 - T1.^2)...
        + coef_CO2(3)./3*(T2.^3 - T1.^3)...
        + coef_CO2(4)./4*(T2.^4 - T1.^4); %[J/mol]
    
    dh_w = coef_w(1).*(T2 - T1) ...
        + coef_w(2)./2*(T2.^2 - T1.^2)...
        + coef_w(3)./3*(T2.^3 - T1.^3)...
        + coef_w(4)./4*(T2.^4 - T1.^4); %[J/mol]
    
    dh_N2 = coef_N2(1).*(T2 - T1) ...
        + coef_N2(2)./2*(T2.^2 - T1.^2)...
        + coef_N2(3)./3*(T2.^3 - T1.^3)...
        + coef_N2(4)./4*(T2.^4 - T1.^4); %[J/mol]
    
    dh_O2 = coef_O2(1).*(T2 - T1) ...
        + coef_O2(2)./2*(T2.^2 - T1.^2)...
        + coef_O2(3)./3*(T2.^3 - T1.^3)...
        + coef_O2(4)./4*(T2.^4 - T1.^4); %[J/mol]
    
    denom = Phi*(H_bar_CO2 + dh_CO2) + c*Phi*(H_bar_H2O + dh_w) ...
        + d* dh_N2 + (a - Phi*(2+c))*dh_O2;
    
    RHS = num/denom;
    
    end
    
    T_a = RHS;
    

end
 