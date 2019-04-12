function [T_1] = varTempEffTurb(T_1s, T_0, eff)
% Outputs the actual or isentropic temperature

T_1 = T_1s;
RHS = dh(T_1s, T_0) ./ dh(T_1, T_0);
dT = 0.001;
for i = 1:2
    while abs(RHS(i) - eff) > 0.001
        if eff > RHS(i)
            T_1(i) = T_1(i) + dT;
        end
        if eff < RHS(i)
            T_1(i) = T_1(i) - dT;
        end
        RHS = dh(T_1s, T_0) ./ dh(T_1, T_0);
    end
end
end