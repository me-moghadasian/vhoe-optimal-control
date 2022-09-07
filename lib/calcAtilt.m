function Atilt = calcAtilt(A, base, n, noT)
    Atilt = [];    
    for k = 1:size(base, 2)
        colind = [];
        for j = 1:size(base, 1)
            if (base(j, k)~=0)
                colind = [colind, j];
            end
        end
        Atilt = [Atilt, makeItDiag(A(:, colind),n , noT)'];
    end
end