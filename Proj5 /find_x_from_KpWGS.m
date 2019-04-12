function x = find_x_from_KpWGS(Kp)
dx = 1e-6;
x = 0;
Kp_est = 0; 

while Kp > Kp_est
    y = x + 3;          % N_H2
    z = 2 - x;          % N_CH4
    w = 1 - x;          % N_H2O

    if (x >= 0 && y >= 0 && w >=0 && z >= 0)
        Kp_est = (x * y)/(w * z);
    else 
        break
    end
    x = x + dx; 
end
end 