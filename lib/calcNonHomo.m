function [nonH, Gamma] = calcNonHomo(PROBLEM, I)
    n = PROBLEM.n;
    c = PROBLEM.c;
    nonH  = sparse(2*n*c, rCombRep(n, I));
    Gamma = [];
    if I==1
        return;
    end
    
    gmaCo = calcGammaCoeff(n, c, I);

    Gamma = cell(I, 1);
    for j = 2:I
        Gamma{j} = 0;
        for k = 1:(I-j+1)
            M   = PROBLEM.Gam{1, k};
            F   = operF(M, PROBLEM.PHIz, j,...
                rCombRep(PROBLEM.n, k), c);
            
            M   = PROBLEM.Gam{j-1, I-k};
            G   = operG(M, PROBLEM.PHIx, I, I-k, c);
            
            Gamma{j} = Gamma{j} + 1.*F*G;
        end
        
        Y     = PROBLEM.D{j}*...
            makeItNonDiag(Gamma{j}, rCombRep(2*PROBLEM.n, j), c);
        nonH  = nonH + Y;
    end
    nonH = nonH.*gmaCo;
end