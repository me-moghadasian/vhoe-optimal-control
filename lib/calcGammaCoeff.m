function gmaCo = calcGammaCoeff(n, c, I)
    W = calcW(n, I);
    Wfac = factorial(W);
    polyCoeff = ones(1, rCombRep(n, I));
    for k = 1:n
        polyCoeff = polyCoeff.*Wfac(k, :);
    end
    polyCoeff = polyCoeff/factorial(I);   
    gmaCo = repmat(polyCoeff, 2*n*c, 1);
end