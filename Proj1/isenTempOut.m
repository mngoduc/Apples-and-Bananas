function [t_out] = isenTempOut(t_in, p_in, p_out)
%Calculates the isentropic temperature out from the input temp, input
%pressure and output pressure based of average k

t_temp = t_in;
for i = 1:10
    [~, k_temp] = specHeatAir(t_temp);
    t_out = t_in .* (p_out ./ p_in) .^ ((k_temp-1)./k_temp); % Isentropic stagnation temperature after fan
    t_temp = t_out;
end
end