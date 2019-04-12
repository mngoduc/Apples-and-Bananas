function T = Tm2T(Tm, RF, M, k)
    T = Tm./(1+RF*((k-1)./2).*M.^2);
end 