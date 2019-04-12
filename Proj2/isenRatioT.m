function ratio = isenRatioT(k, M)
    % using isentropic ratios, returns To/T
    ratio = 1 + (k-1)./2.*M.^2;
end