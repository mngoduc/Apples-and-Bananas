function ratio = isenRatioP(k, M)
    % using isentropic ratios, returns Po/P
    ratio = (1 + (k-1)./2.*M.^2).^(k./(k-1));
end