function res = makeItDiag(A, noVar, m)
    colNo = size(A, 2); 
    res = zeros(noVar*m, colNo*m);
    for i1 = 1:noVar
        for i2 = 1:colNo
            res([1:m]+(i1-1)*m, [1:m]+(i2-1)*m) = ...
                diag(A([1:m]+(i1-1)*m, i2));
        end
    end
end