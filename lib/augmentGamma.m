function Gamma = augmentGamma(PROBLEM, k, Gam)
    gmaCo = calcGammaCoeff(PROBLEM.n, PROBLEM.c, k);

    Gamma = PROBLEM.Gam;
    Gamma{1, k} = sparse(makeItDiag(PROBLEM.S{k}./gmaCo,...
        2*PROBLEM.n, PROBLEM.c));
    for j = 2:k
        Gamma{j, k} = sparse(Gam{j});
    end
end