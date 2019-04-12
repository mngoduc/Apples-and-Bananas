function [MFP,Po,M] = ToP2MFP(mDot,R,To,A,k,P)
    % ONLY TAKES SCALARS
    M = 0; % initilize mach
    dM = 0.00001;
    Po = P;    % initial guess for Po
    MFP = mDot*sqrt(R*To)/(A*Po);
    LHS = MFP;
    RHS = M*sqrt(k)*(1+(k-1)/2*M^2)^-((k+1)/(2*(k-1)));
    while LHS > RHS % iterate until found Po and M that satisfy MFP eqn
        M = M + dM;
        Po = P*isenRatioP(k, M);
        LHS = mDot*sqrt(R*To)/(A*Po);
        RHS = M*sqrt(k)*(1+(k-1)/2*M^2)^-((k+1)/(2*(k-1)));
    end
    MFP = LHS;
end