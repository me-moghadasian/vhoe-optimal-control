function res = makeItNonDiag(A, noVar, noT)
    colNo = size(A, 2)/noT;
    res = zeros(noVar*noT, colNo);
    for i1 = 1:noVar
        for i2 = 1:colNo
            res((i1-1)*noT+1:i1*noT, i2) = ...
                diag(A((i1-1)*noT+1:i1*noT, (i2-1)*noT+1:i2*noT));
        end
    end
end