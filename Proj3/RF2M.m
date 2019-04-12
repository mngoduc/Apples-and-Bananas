function M = RF2M(Tm,To,RF,k)
    % ONLY ACCEPTS SCALARS
    M = 0;
    dm = 0.001;
    LHS = Tm/To;
    RHS = (1+RF*((k-1)/2)*M^2)/(1+((k-1)/2)*M^2);
    while LHS > RHS
        M = M + dm;
        RHS = (1+RF*((k-1)/2)*M^2)/(1+((k-1)/2)*M^2);
    end
end