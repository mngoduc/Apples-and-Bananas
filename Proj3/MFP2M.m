function  M = MFP2M(k,MFP)
    % using MFP equation, iterates to find mach number
    % ONLY ACCEPTS SCALARS
    M = fzero(@(M) M*sqrt(k)*(1+(k-1)/2*M^2)^-((k+1)/(2*(k-1)))-MFP, 0.25);
end