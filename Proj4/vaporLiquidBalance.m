function [beta, gamma] = vaporLiquidBalance (T, y_test, y_max, N_a, N_H2O, alpha)

    if (length(y_test) == 1)
        for i = 1:max(length(T),length(y_max))
            if y_test > y_max(i)
                y_actual(i) = y_max(i);
                N_v(i) = y_max(i)*N_a/(1-y_max(i));
                N_l(i) = N_H2O - N_v(i); 
            else % y_test <= y_max
                y_actual(i) = y_test;
                N_v(i) = N_H2O ;
                N_l(i) = 0; 
            end
        end
    elseif (length(y_max) == 1)    
        for i = 1:length(y_test)
            if y_test(i) > y_max
                y_actual(i) = P2.y_max;
                N_v(i) = y_max * N_a./(1-y_max);
                N_l(i) = N_H2O - N_v(i); 
            else % y_test(i) <= y_max
                y_actual(i) = y_test(i);
                N_v(i) = N_H2O ;
                N_l(i) = 0; 
            end
        end
    else
        for i = 1:length(y_test)
            if y_test(i) > y_max(i)
                y_actual(i) = y_max(i);
                N_v(i) = y_max(i) * N_a/(1-y_max(i));
                N_l(i) = N_H2O(i) - N_v(i); 
            else % y_test <= y_max
                y_actual(i) = y_test(i);
                N_v(i) = N_H2O(i);
                N_l(i) = 0; 
            end
        end
    end
% beta = moles of water vapor 
% gamma = moles of liquid water
beta = N_v; 
gamma = N_l; 
end